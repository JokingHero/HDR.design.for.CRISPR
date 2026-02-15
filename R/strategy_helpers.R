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

    # Check Genomic Overlap
    if (length(findOverlaps(candidate_mut, selected_gr)) > 0) {
      next
    }

    # Check Codon Overlap
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
#' @description Pre-calculates standardized metrics on the var_data GRanges to
#'   facilitate downstream sorting. Handles missing data (e.g., non-model organisms)
#'   by imputing neutral defaults.
#' @param var_data The main GRanges object of candidate SNPs.
#' @param optimization_scheme Typically "balanced"
#' @param benign_cadd_threshold Numeric. CADD score below which a variant is
#'   considered high-confidence benign. Default 15.
#' @param ag_threshold Numeric. AlphaGenome Composite Score above which a variant
#'   is considered risky. Default 1.5 and implies strong signal in at least one
#'   modality and moderate support in others. Research suggests >1.0 is moderate risk.
#' @return The augmented var_data GRanges object with standardized columns:
#'   `penalty_score`, `disruption_tier`, `safety_tier`, etc.
#' @keywords internal
#'
augment_var_data_with_scores <- function(var_data,
                                         optimization_scheme = "balanced",
                                         benign_cadd_threshold = 15,
                                         ag_threshold = 1.5) {
  if (length(var_data) == 0) return(var_data)

  # ============================================================================
  # 1. GUIDE POSITION (Guaranteed Feature)
  # ============================================================================
  var_data$disruption_tier <- 5
  var_data$disruption_tier[var_data$position_in_guide >= 10] <- 4
  var_data$disruption_tier[var_data$position_in_guide >= 13] <- 3
  var_data$disruption_tier[var_data$position_in_guide >= 17] <- 2
  var_data$disruption_tier[var_data$position_in_guide >= 22] <- 1

  # ============================================================================
  # 2. DBSNP / KNOWN VARIATION (Optional Feature)
  # ============================================================================
  var_data$dbSNP_priority <- 3

  if (!is.null(var_data$dbSNP)) {
    is_known_variant <- sapply(var_data$dbSNP, function(df) {
      if (nrow(df) > 0) any(df$is_known_variant) else NA
    })
    var_data$dbSNP_priority[!is.na(is_known_variant) & !is_known_variant] <- 2
    var_data$dbSNP_priority[is_known_variant] <- 1
  }

  # ============================================================================
  # 3. NON-CODING OVERLAPS (Guaranteed Annotation Feature)
  # ============================================================================
  var_data$has_nc_overlap <- FALSE
  if (!is.null(var_data$noncoding)) {
    var_data$has_nc_overlap <- sapply(var_data$noncoding, nrow) > 0
  }

  # ============================================================================
  # 4. CADD SCORES (Optional Feature)
  # ============================================================================
  if (is.null(var_data$CADD)) var_data$CADD <- NA
  var_data$cadd_imputed <- var_data$CADD
  var_data$cadd_imputed[is.na(var_data$cadd_imputed)] <- 15

  # ============================================================================
  # 5. ALPHAGENOME (Splicing Composite Score)
  # ============================================================================
  # Research Update: We now use the Composite Score (Sites + Usage + Junctions/5).
  # This single number captures magnitude, motif disruption, and isoform switching.

  ag_score <- rep(0, length(var_data))
  ag_has_data <- rep(FALSE, length(var_data))

  if ("alphagenome_composite_score" %in% names(mcols(var_data))) {
    raw <- var_data$alphagenome_composite_score
    ag_has_data <- !is.na(raw)
    ag_score[ag_has_data] <- raw[ag_has_data]
  } else if ("ag_composite_score" %in% names(mcols(var_data))) {
    raw <- var_data$ag_composite_score
    ag_has_data <- !is.na(raw)
    ag_score[ag_has_data] <- raw[ag_has_data]
  } else {
    # Fallback for legacy data or if composite calculation failed:
    # Try to sum individually if columns exist
    val_sites <- if ("alphagenome_sites_max" %in% names(mcols(var_data))) {
      var_data$alphagenome_sites_max
    } else if ("alphagenome_SPLICE_SITES" %in% names(mcols(var_data))) {
      var_data$alphagenome_SPLICE_SITES
    } else rep(NA_real_, length(var_data))

    val_usage <- if ("alphagenome_usage_max" %in% names(mcols(var_data))) {
      var_data$alphagenome_usage_max
    } else if ("alphagenome_SPLICE_SITE_USAGE" %in% names(mcols(var_data))) {
      var_data$alphagenome_SPLICE_SITE_USAGE
    } else rep(NA_real_, length(var_data))

    val_junctions <- if ("alphagenome_junctions_max" %in% names(mcols(var_data))) {
      var_data$alphagenome_junctions_max
    } else if ("alphagenome_SPLICE_JUNCTIONS" %in% names(mcols(var_data))) {
      var_data$alphagenome_SPLICE_JUNCTIONS
    } else rep(NA_real_, length(var_data))

    ag_has_data <- !(is.na(val_sites) & is.na(val_usage) & is.na(val_junctions))

    # Note: treating missing components as 0 only when at least one AG component exists.
    val_sites[is.na(val_sites)] <- 0
    val_usage[is.na(val_usage)] <- 0
    val_junctions[is.na(val_junctions)] <- 0
    ag_score <- val_sites + val_usage + (val_junctions / 5)
  }
  var_data$ag_composite_score <- ag_score
  var_data$ag_impact_score <- ag_score
  var_data$ag_has_data <- ag_has_data

  # ============================================================================
  # 6. FINAL SAFETY TIER (The "Meta" Feature)
  # ============================================================================
  # 1 = Proven Benign (dbSNP Match)
  # 2 = High Confidence Benign (Low CADD + Low AG + No NC Overlap)
  # 3 = Neutral / No Data (Imputed CADD + No AG + No NC Overlap)
  # 4 = Predicted Risky (High CADD or High AG)
  # 5 = Structural/Annotation Risk (Non-coding overlap)
  tier <- rep(3, length(var_data))

  # Risk Factors
  is_ag_risky <- ag_score > ag_threshold
  is_risky_pred <- (var_data$cadd_imputed > benign_cadd_threshold) | is_ag_risky
  tier[is_risky_pred] <- 4
  tier[var_data$has_nc_overlap] <- 5

  # Safety Factors
  has_real_data <- !is.na(var_data$CADD)
  is_low_cadd <- var_data$cadd_imputed < benign_cadd_threshold

  # Predicted Safe: Needs Low CADD AND Low AG AND No Overlap
  tier[has_real_data & is_low_cadd & !is_risky_pred & !var_data$has_nc_overlap] <- 2

  # Known Benign Overrides predictions
  tier[var_data$dbSNP_priority == 1 & !var_data$has_nc_overlap] <- 1
  var_data$safety_tier <- tier

  # ============================================================================
  # 7. CALCULATE PRIORITY GROUP (Switch Case)
  # ============================================================================

  prio <- rep(99, length(var_data))
  d_tier <- var_data$disruption_tier
  s_tier <- var_data$safety_tier

  if (optimization_scheme == "balanced") {
    # 1. Platinum: PAM/Crit (1-2) + Benign/HighConf (1-2)
    prio[d_tier <= 2 & s_tier <= 2] <- 1
    # 2. Gold: Seed (3) + Benign/HighConf (1-2)
    prio[d_tier == 3 & s_tier <= 2] <- 2
    # 3. Silver: PAM/Crit (1-2) + Neutral (3)
    prio[d_tier <= 2 & s_tier == 3] <- 3
    # 4. Bronze: Seed (3) + Neutral (3)
    prio[d_tier == 3 & s_tier == 3] <- 4
    # 5. Iron: Distal (4-5) + Benign/HighConf (1-2)
    prio[d_tier >= 4 & s_tier <= 2] <- 5
    # 6. Lead: Distal (4-5) + Neutral (3)
    prio[d_tier >= 4 & s_tier == 3] <- 6
    # 7. Maverick: PAM (1-2) + Risky (4)
    prio[d_tier <= 2 & s_tier == 4] <- 7
    # 8. Radioactive: Everything else
    prio[s_tier >= 4 & d_tier >= 3] <- 8
    prio[s_tier == 5] <- 9

  } else if (optimization_scheme == "disruption_first") {
    # 1. Guaranteed Kill: PAM (1-2)
    prio[d_tier <= 2 & s_tier <= 4] <- 1
    # 2. Likely Kill: Seed (3)
    prio[d_tier == 3 & s_tier <= 4] <- 2
    # 3. Weak: Distal (4-5)
    prio[d_tier >= 4 & s_tier <= 4] <- 3
    # 4. Invalid
    prio[s_tier == 5] <- 4
  } else if (optimization_scheme == "safety_first") {
    prio[s_tier == 1] <- 1
    prio[s_tier == 2] <- 2
    prio[s_tier == 3] <- 3
    prio[s_tier == 4] <- 4
    prio[s_tier == 5] <- 5
  } else {
    stop("Invalid optimization_scheme.")
  }

  var_data$priority_group <- prio
  return(var_data)
}

#' @title Find the best set of SNPs for a single guide
#' @description Sorts candidate SNPs based on a chosen scheme and greedily
#'   selects the top N non-overlapping ones.
#' @return A GRanges object of the selected mutations.
#' @keywords internal
#'
find_best_snps_for_guide <- function(guide_snps, N, optimization_scheme) {
  if (length(guide_snps) == 0) return(GRanges())

  ordering <- switch(
    optimization_scheme,
    "balanced" = order(
      guide_snps$priority_group,
      guide_snps$disruption_tier,
      guide_snps$cadd_imputed,
      guide_snps$ag_impact_score,
      -guide_snps$position_in_guide
    ),
    "disruption_first" = order(
      guide_snps$priority_group,
      guide_snps$safety_tier,
      guide_snps$cadd_imputed,
      guide_snps$ag_impact_score
    ),
    "safety_first" = order(
      guide_snps$priority_group,
      guide_snps$disruption_tier,
      guide_snps$ag_impact_score,
      -guide_snps$position_in_guide
    ),
    stop("Invalid optimization_scheme")
  )

  sorted_snps <- guide_snps[ordering]
  select_non_overlapping_mutations(sorted_snps, N)
}

#' @title Assess Guide Disruption on Final Template
#' @description Aligns the guide sequence (with PAM) to the final generated template
#' sequence to precisely calculate disruption metrics.
#' @keywords internal
#'
assess_guide_disruption <- function(guide, template_seq) {
  aln <- pwalign::pairwiseAlignment(
    DNAString(guide$with_pam), template_seq,
    type = "global-local",
    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
      match = 1, mismatch = 0, baseOnly = F),
    gapOpening = 0, gapExtension = -1
  )

  # Get the aligned pattern (guide) and subject (template) substrings
  # Note: pattern(aln) returns the aligned version (with gaps), but we want the substring
  # corresponding to the alignment on the subject.

  is_plus <- as.character(guide$strand) == "+"

  # Extract aligned strings (these are characters)
  # We need to look at mismatches.
  # `pattern(aln)` is the guide (query). `subject(aln)` is the template part.
  # `compareStrings` gives a string with "?" for mismatch.
  comp <- pwalign::compareStrings(aln)

  # Parse mismatches
  mismatches <- unlist(strsplit(comp, ""))
  is_mismatch <- mismatches %in% c("?", "+") # + is gap

  # Mismatch indices (1-based relative to the 23bp guide fragment)
  mm_indices <- which(is_mismatch)
  list(
    total_disrupton = length(mm_indices),
    pam_disrupted = sum(mm_indices %in% if (is_plus) c(22, 23) else c(1, 2)),
    seed_disrupted = sum(mm_indices %in% if (is_plus) c(11:20) else c(4:13)),
    aln_guide = as.character(pwalign::pattern(aln)),
    aln_template = as.character(pwalign::subject(aln))
  )
}

#' @title Create HDR template and probes from a set of mutations
#' @description Injects selected SNPs into the base template sequence, calculates
#'   summary statistics based on re-alignment, and designs validation probes.
#' @return A list containing the 'template' and 'probes' GRanges.
#' @keywords internal
#'
create_template_and_probes <- function(selected_muts,
                                       variants_in_editw,
                                       variants_genomic_on_ts,
                                       design_id,
                                       edit_region,
                                       source_genomic_seq,
                                       current_guide,
                                       do_probes,
                                       probe_params) {

  muts_to_inject <- variants_in_editw
  snvs_introduced <- ""
  total_cadd <- 0
  max_ag_score <- NA_real_
  total_snp_quality <- 0
  any_overlaps_nc <- FALSE

  if (length(selected_muts) > 0) {
    # Map selected SNPs to edit window
    muts_in_editw <- pmapToTranscripts(selected_muts, edit_region)
    muts_in_editw$ALT <- selected_muts$ALT

    # Combine introduced SNPs with original template variants for injection
    muts_to_inject <- sort(c(variants_in_editw, muts_in_editw))

    # Stats
    snvs_introduced <- paste0(names(selected_muts), collapse = ";")
    total_cadd <- sum(selected_muts$cadd_imputed, na.rm = TRUE)
    total_snp_quality <- sum(selected_muts$dbSNP_priority, na.rm = TRUE)
    any_overlaps_nc <- any(selected_muts$has_nc_overlap)
    if ("ag_composite_score" %in% names(mcols(selected_muts))) {
      has_ag_data <- if ("ag_has_data" %in% names(mcols(selected_muts))) {
        !is.na(selected_muts$ag_has_data) & selected_muts$ag_has_data
      } else {
        !is.na(selected_muts$ag_composite_score)
      }

      if (any(has_ag_data)) {
        max_ag_score <- max(
          selected_muts$ag_composite_score[has_ag_data], na.rm = TRUE)
        if (!is.finite(max_ag_score)) max_ag_score <- NA_real_
      }
    } else if ("alphagenome_composite_score" %in% names(mcols(selected_muts))) {
      has_ag_data <- !is.na(selected_muts$alphagenome_composite_score)
      if (any(has_ag_data)) {
        max_ag_score <- max(
          selected_muts$alphagenome_composite_score[has_ag_data], na.rm = TRUE)
        if (!is.finite(max_ag_score)) max_ag_score <- NA_real_
      }
    }
  }

  # Inject all mutations to create the final HDR sequence
  final_hdr_seq <- if (length(muts_to_inject) > 0) {
    replaceAt(source_genomic_seq,
              at = ranges(muts_to_inject),
              value = DNAStringSet(muts_to_inject$ALT))
  } else {
    source_genomic_seq
  }

  # 3. Score Guide against Final Template
  # This replaces the old static calculation
  guide_stats <- assess_guide_disruption(current_guide, final_hdr_seq)

  # 4. Create Template Object
  hdr_template_coord_map <- build_variant_layout(
    muts_to_inject, nchar(source_genomic_seq))
  template_coords <- remap_target_to_genomic(
    target = GRanges("target_seq", ranges = 1:width(edit_region)),
    coordinate_map = hdr_template_coord_map,
    window_genomic = edit_region,
    variants_genomic = c(selected_muts, variants_genomic_on_ts))

  template_gr <- edit_region
  names(template_gr) <- design_id
  template_gr$coords <- template_coords$coords
  template_gr$sequence <- as.character(final_hdr_seq)
  template_gr$snvs_introduced <- snvs_introduced
  template_gr$pam_disrupted_count <- guide_stats$pam_disrupted
  template_gr$seed_disrupted_count <- guide_stats$seed_disrupted
  template_gr$total_disruption_count <- guide_stats$total_disrupton
  template_gr$aln_guide <- guide_stats$aln_guide
  template_gr$aln_template <- guide_stats$aln_template
  template_gr$total_cadd <- total_cadd
  template_gr$max_alphagenome_score <- max_ag_score
  template_gr$any_overlaps_noncoding <- any_overlaps_nc
  template_gr$total_snp_quality_score <- total_snp_quality

  # --- 5. Design Probes if requested ---
  probes_out <- GRanges()
  # If there are no muts to cover, no point in HDR probes
  if (do_probes) {
    candidates <- design_probes(
      s = final_hdr_seq,
      genomic_context = edit_region,
      coordinate_map = hdr_template_coord_map,
      variants_genomic = c(selected_muts, variants_genomic_on_ts),
      tmin = probe_params$tmin, tmax = probe_params$tmax,
      len_min = probe_params$len_min, len_max = probe_params$len_max)

    if (nrow(candidates) > 0) {
      mcols(selected_muts) <- NULL
      mcols(variants_genomic_on_ts) <- NULL
      probes_out <- select_probes(
        muts_to_cover = c(selected_muts, variants_genomic_on_ts),
        candidates = candidates,
        temp_name = design_id)

      if (nrow(probes_out) > 0) {
        probes_out$names <- paste0("HDR_probe_", seq_len(nrow(probes_out)), "_for_", design_id)
      }
    }
  }

  list(template = template_gr, probes = probes_out)
}
