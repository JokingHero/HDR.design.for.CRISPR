#' Selects mutations that are non-overlapping at both the genomic and codon level.
#'
#' @param gr A GRanges object of candidate mutations, assumed to be pre-sorted by priority.
#' @param N The maximum number of mutations to select.
#' @return A GRanges object containing the top N non-overlapping mutations.
#' @keywords internal
#'
select_non_overlapping_mutations <- function(gr, N) {
  if (length(gr) == 0) {
    return(gr[0])
  }

  selected_gr <- gr[0]
  seen_codon_keys <- c()
  for (i in seq_along(gr)) {
    if (length(selected_gr) >= N) {
      break
    }

    candidate_mut <- gr[i]
    if (length(findOverlaps(candidate_mut, selected_gr)) > 0) {
      next
    }

    current_cds_df <- candidate_mut$CDS[[1]]
    current_keys <- c()
    has_codon_data <- !is.null(current_cds_df) && nrow(current_cds_df) > 0

    if (has_codon_data) {
      current_keys <- paste(current_cds_df$tx_id, current_cds_df$codon_num, sep = "_")
      if (any(current_keys %in% seen_codon_keys)) {
        next
      }
    }

    selected_gr <- c(selected_gr, candidate_mut)
    if (has_codon_data) {
      seen_codon_keys <- union(seen_codon_keys, current_keys)
    }
  }

  return(selected_gr)
}

#' @title Augment var_data with scoring and priority columns
#' @description Pre-calculates several columns on the var_data GRanges to
#'   facilitate the different optimization sorting schemes.
#' @param var_data The main GRanges object of candidate SNPs.
#' @param pam_proximal_threshold An integer defining what `position_in_guide`
#'   is considered PAM-proximal.
#' @param benign_cadd_threshold A numeric CADD score below which a variant is
#'   considered "predicted benign".
#' @return The augmented var_data GRanges object.
#' @keywords internal
augment_var_data_with_scores <- function(var_data,
                                         pam_proximal_threshold = 12,
                                         benign_cadd_threshold = 5) {
  if (length(var_data) == 0) return(var_data)

  # --- 1. Standardize Safety/Confidence Metrics ---
  # dbSNP Priority (lower is better)
  var_data$dbSNP_priority <- 3 # Default: Unknown
  if (!is.null(var_data$dbSNP)) {
    is_known_variant <- sapply(var_data$dbSNP, function(df) {
      if (nrow(df) > 0) any(df$is_known_variant) else NA
    })
    var_data$dbSNP_priority[!is.na(is_known_variant) & !is_known_variant] <- 2 # Known, but not matching allele
    var_data$dbSNP_priority[is_known_variant] <- 1 # Gold standard: Known benign
  }

  # Non-coding Overlap
  var_data$has_nc_overlap <- FALSE
  if (!is.null(var_data$noncoding)) {
    var_data$has_nc_overlap <- sapply(var_data$noncoding, nrow) > 0
  }

  # Ensure CADD exists
  if (is.null(var_data$CADD)) var_data$CADD <- 0

  # --- 2. Calculate Tiers for "Balanced" Scheme ---
  var_data$is_pam_proximal <- var_data$position_in_guide > pam_proximal_threshold
  is_known_benign <- var_data$dbSNP_priority == 1
  is_predicted_benign <- var_data$CADD < benign_cadd_threshold & !var_data$has_nc_overlap

  # Assign tiers from worst to best (lower tier number is better)
  var_data$priority_tier <- 6 # Tier 6: Fallback (default)
  var_data$priority_tier[var_data$is_pam_proximal & !is_known_benign & !is_predicted_benign] <- 5
  var_data$priority_tier[is_predicted_benign & !var_data$is_pam_proximal & !is_known_benign] <- 4
  var_data$priority_tier[is_predicted_benign & var_data$is_pam_proximal & !is_known_benign] <- 3
  var_data$priority_tier[is_known_benign & !var_data$is_pam_proximal] <- 2
  var_data$priority_tier[is_known_benign & var_data$is_pam_proximal] <- 1 # Tier 1: Perfect
  return(var_data)
}

#' @title Find the best set of SNPs for a single guide
#' @description Sorts candidate SNPs based on a chosen scheme and greedily
#'   selects the top N non-overlapping ones.
#' @return A GRanges object of the selected mutations.
#' @keywords internal
find_best_snps_for_guide <- function(guide_snps, N, optimization_scheme) {
  if (length(guide_snps) == 0) return(GRanges())

  ordering <- switch(
    optimization_scheme,
    "balanced" = order(
      guide_snps$priority_tier,      # Tier 1 is best
      -guide_snps$position_in_guide, # Within tier, proximal is best
      guide_snps$CADD                # Further tie-breaker
    ),
    "safety_first" = order(
      guide_snps$dbSNP_priority,       # Known benign is best
      guide_snps$has_nc_overlap,       # No overlap is best
      guide_snps$CADD,                 # Low CADD is best
      -guide_snps$position_in_guide    # Tie-breaker: disruption
    ),
    "disruption_first" = order(
      -guide_snps$position_in_guide,   # Proximal is best
      guide_snps$dbSNP_priority,       # Tie-breaker: safety
      guide_snps$has_nc_overlap,
      guide_snps$CADD
    ),
    stop("Invalid optimization_scheme")
  )

  sorted_snps <- guide_snps[ordering]
  select_non_overlapping_mutations(sorted_snps, N)
}

#' @title Find all valid combinations of SNPs for a single guide
#' @description Generates all combinations of a given size and filters them
#'   for validity (no overlaps).
#' @return A list of GRanges objects, each a valid combination.
#' @keywords internal
find_all_snp_combinations_for_guide <- function(guide_snps, mpt) {
  if (length(guide_snps) < mpt) return(list())

  combs_indices <- combn(seq_along(guide_snps), mpt, simplify = FALSE)

  valid_combinations <- lapply(combs_indices, function(indices) {
    candidate_muts <- guide_snps[indices]
    # Use select_non_overlapping_mutations as a validity check
    selected <- select_non_overlapping_mutations(candidate_muts, mpt)
    if (length(selected) == mpt) {
      return(selected)
    } else {
      return(NULL)
    }
  })

  # Remove NULLs from the list
  valid_combinations[!sapply(valid_combinations, is.null)]
}

#' @title Create HDR template and probes from a set of mutations
#' @description Injects selected SNPs into the base template sequence, calculates
#'   summary statistics, and designs validation probes.
#' @return A list containing the 'template' and 'probes' GRanges.
#' @keywords internal
#'
create_template_and_probes <- function(selected_muts,
                                       variants_in_editw,
                                       variants_genomic_on_ts,
                                       design_id,
                                       edit_region,
                                       source_genomic_seq,
                                       do_probes,
                                       probe_params) {
  # Return empty list if no mutations were selected
  if (is.null(selected_muts) || length(selected_muts) == 0) {
    return(list(template = GRanges(), probes = GRanges()))
  }

  # Combine introduced SNPs with original template variants for injection
  muts_in_editw <- pmapToTranscripts(selected_muts, edit_region)
  muts_in_editw$ALT <- selected_muts$ALT
  all_muts_in_editw <- sort(c(variants_in_editw, muts_in_editw))

  # Inject all mutations to create the final HDR sequence
  final_hdr_seq <- replaceAt(source_genomic_seq,
                             at = ranges(all_muts_in_editw),
                             value = DNAStringSet(all_muts_in_editw$ALT))

  # --- Summary Statistics ---
  total_cadd <- sum(selected_muts$CADD, na.rm = TRUE)
  pam_disrupted_count <- sum(selected_muts$position_in_guide >= 22)
  guide_disrupted_count <- length(unique(selected_muts$guide_name))
  any_overlaps_noncoding <- any(selected_muts$has_nc_overlap)
  total_snp_quality_score <- sum(selected_muts$dbSNP_priority)

  # --- Create the final template GRanges object ---
  template_gr <- edit_region
  names(template_gr) <- design_id
  template_gr$sequence <- as.character(final_hdr_seq)
  template_gr$snps_introduced <- paste(names(selected_muts), collapse = ";")
  template_gr$pam_disrupted_count <- pam_disrupted_count
  template_gr$guide_disrupted_count <- guide_disrupted_count
  template_gr$total_cadd <- total_cadd
  template_gr$any_overlaps_noncoding <- any_overlaps_noncoding
  template_gr$total_snp_quality_score <- total_snp_quality_score

  # --- Design Probes if requested ---
  probes_out <- GRanges()
  if (do_probes) {
    # We need a coordinate map for the final_hdr_seq
    hdr_template_coord_map <- build_variant_layout(
      all_muts_in_editw, nchar(source_genomic_seq))
    all_genomic_variants <- c(selected_muts, variants_genomic_on_ts)

    candidates <- design_probes(
      s = final_hdr_seq,
      genomic_context = edit_region,
      coordinate_map = hdr_template_coord_map,
      variants_genomic = all_genomic_variants,
      tmin = probe_params$tmin, tmax = probe_params$tmax,
      len_min = probe_params$len_min, len_max = probe_params$len_max)

    if (nrow(candidates) > 0) {
      probes_out <- select_probes(
        muts_to_cover = selected_muts,
        candidates = candidates,
        temp_name = design_id)

      if (nrow(probes_out) > 0) {
        probes_out$names <- paste0("HDR_probe_", seq_len(nrow(probes_out)), "_for_", design_id)
      }
    }
  }

  list(template = template_gr, probes = probes_out)
}
