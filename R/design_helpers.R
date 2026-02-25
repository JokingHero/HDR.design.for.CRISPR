#' Primer 3
#' @description Design primers with the use of primer 3
#' We want at least one primer to be outside of the template_range
#' We use `template_range_extended` as the region where we can design
#' We check for off-targets
#' We want to make sure the target_on_template_extended is captured inside the PCRed sequence
#' @keywords internal
design_primers <- function(
    template_range_extended, template_range, variant_with_allowed, genome, primer3_path){
  message("Designing primers...")
  options <- c(
    "PRIMER_SEQUENCE_ID",
    "SEQUENCE_TEMPLATE",
    "SEQUENCE_TARGET",
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST",
    "PRIMER_PICK_LEFT_PRIMER",
    "PRIMER_PICK_INTERNAL_OLIGO",
    "PRIMER_PICK_RIGHT_PRIMER",
    "PRIMER_OPT_SIZE",
    "PRIMER_MIN_SIZE",
    "PRIMER_MAX_SIZE",
    "PRIMER_MAX_NS_ACCEPTED",
    "PRIMER_PRODUCT_SIZE_RANGE",
    "P3_FILE_FLAG",
    "PRIMER_EXPLAIN_FLAG")
  template_range_extended_seq <- getSeq(genome, template_range_extended)[[1]]
  variant_with_allowed_template_extended <- pmapToTranscripts(
    variant_with_allowed, template_range_extended)
  template_range_on_template_extended <- pmapToTranscripts(
    template_range, template_range_extended)
  values <- c(
    "template",
    as.character(template_range_extended_seq),
    paste0(start(variant_with_allowed_template_extended), ",",
           width(variant_with_allowed_template_extended)),
    paste0(1, ",", start(template_range_on_template_extended), ",, ; ,,",
           end(template_range_on_template_extended), ",",
           nchar(template_range_extended_seq) - end(template_range_on_template_extended) - 1),
    "1", "0", "1",
    getOption("PRIMER_OPT_SIZE", default = "22"),
    getOption("PRIMER_MIN_SIZE", default = "18"),
    getOption("PRIMER_MAX_SIZE", default = "25"),
    getOption("PRIMER_MAX_NS_ACCEPTED", default = "0"),
    paste0(getOption("PRODUCT_SIZE_MIN", default = "100"), "-",
           getOption("PRODUCT_SIZE_MAX", default = "250")),
    "0", "1"
  )
  opts <- c(paste0(options, "=", values), "=")

  out <- tryCatch({
    system2(command = primer3_path, input = opts, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    warning("Failed to execute primer3. Please check if 'primer3_core' is installed and the path is correct. Error: ", e$message)
    return(GRanges())
  })

  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    warning("Primer3 execution failed with exit code ", status, ". Output:\n", paste(out, collapse="\n"))
    return(GRanges())
  }

  if (is.null(out) || length(out) < 2) {
    warning("Primer3 did not return valid output.")
    return(GRanges())
  }

  out <- strsplit(out, "=")
  out <- out[1:(length(out) - 1)]
  tags <- sapply(out, `[[`, 1)
  res <- sapply(out, `[[`, 2)

  pair_count <- as.numeric(res[tags == "PRIMER_PAIR_NUM_RETURNED"])
  if (is.na(pair_count) || pair_count == 0) {
    message("Primer3 did not find any suitable primer pairs.")
    return(GRanges())
  }
  primers <- list()
  for (i in seq_len(pair_count)) {
    sl <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_SEQUENCE")]
    sr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_SEQUENCE")]
    ll <- res[tags == paste0("PRIMER_LEFT_", i - 1)] #  "start,width"
    lr <- res[tags == paste0("PRIMER_RIGHT_", i - 1)]
    ll <- strsplit(ll, ",")[[1]]
    lr <- strsplit(lr, ",")[[1]]
    lw <- as.numeric(ll[2])
    rw <- as.numeric(lr[2])
    ll <- as.numeric(ll[1]) # coordinates already in BED format with -1
    lr <- as.numeric(lr[1])
    lr <- lr - rw + 1
    # left primer
    # template_range_extended_seq[(ll + 1):(ll + lw)]
    # right primer
    # template_range_extended_seq[(lr + 1):(lr + rw)]

    tmr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_TM")]
    tml <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_TM")]
    gcl <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_GC_PERCENT")]
    gcr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_GC_PERCENT")]
    ps <- res[tags == paste0("PRIMER_PAIR_", i - 1, "_PRODUCT_SIZE")]

    # offtargets detection with up to 2 mismatch
    ot_count <- 0
    for (ch in seqlevels(genome)) {
      # plus strand
      ot <- Biostrings::matchLRPatterns(
        sl, as.character(reverseComplement(DNAString(sr))),
        max.Lmismatch = 2, max.Rmismatch = 2,
        max.gaplength = 500, subject = genome[[ch]])
      # we need to filter too small products
      # we need to filter out repeated overlaps
      ot <- reduce(ranges(ot[nchar(ot) > (lw + rw)]))
      if (as.vector(seqnames(template_range_extended)) == ch) {
        ot <- setdiff(ot, ranges(template_range_extended))
      }
      ot_count <- ot_count + length(ot)

      ot <- Biostrings::matchLRPatterns(
        sr, as.character(reverseComplement(DNAString(sl))),
        max.Lmismatch = 2, max.Rmismatch = 2,
        max.gaplength = 500, subject = genome[[ch]])
      ot <- reduce(ranges(ot[nchar(ot) > (lw + rw)]))
      if (as.vector(seqnames(template_range_extended)) == ch) {
        ot <- setdiff(ot, ranges(template_range_extended))
      }
      if (length(ot) > 0) message(ch, "-")
      ot_count <- ot_count + length(ot)
    }

    primers[[i]] <- c(ps, ot_count, sl, ll, lw, tml, gcl, sr, lr, rw, tmr, gcr)
  }
  primers <- as.data.frame(do.call(rbind, primers))
  colnames(primers) <- c(
    "product_size", "offtarget_count",
    "LEFT_sequence", "LEFT_start_pos", "LEFT_size",
    "LEFT_melting_temp", "LEFT_GC",
    "RIGHT_sequence", "RIGHT_start_pos", "RIGHT_size",
    "RIGHT_melting_temp", "RIGHT_GC")
  primers <- primers[order(primers$offtarget_count, primers$product_size), ]
  primers_genomic_l <- pmapFromTranscripts(IRanges(start = as.numeric(primers$LEFT_start_pos) + 1,
                                                   width = as.numeric(primers$LEFT_size)),
                                           template_range_extended)
  primers_genomic_r <- pmapFromTranscripts(IRanges(start = as.numeric(primers$RIGHT_start_pos) + 1,
                                                   width = as.numeric(primers$RIGHT_size)),
                                           template_range_extended)
  primers$LEFT_genomic_start <- start(primers_genomic_l)
  primers$RIGHT_genomic_start <- start(primers_genomic_r)
  primers$chrom <- as.character(seqnames(primers_genomic_l))
  primers$rownames <- paste0("PAIR ", seq_along(primers_genomic_l))
  return(primers)
}

#' @title Helper write components
#' @description Helper to write CSV and GFF3 for a given component
#' @keywords internal
#'
write_component_files <- function(
    output_dir, design_name, component, component_name) {
  if (methods::is(component, "GRanges")) {
    component <- as.data.frame(component)
    component$coords <- paste0(
      component$seqnames, ":", component$start, "-", component$end)
  }

  if (nrow(component) == 0) {
    return(invisible(NULL))
  }

  component$seqnames <- NULL
  component$start <- NULL
  component$end <- NULL
  component$width <- NULL
  component$strand <- NULL
  component$names <- rownames(component)
  front <- intersect(c("names", "coords"), names(component))
  component <- component[, c(front, setdiff(names(component), front)), drop = FALSE]

  # Write 1-based coordinate CSV file
  write.table(
    component,
    file = file.path(output_dir, paste0(design_name, "_", component_name, "_1based.csv")),
    quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
}

#' Calculate Pareto Rank (Skyline)
#'
#' Performs non-dominated sorting (Pareto ranking) on a set of candidates.
#' Rank 1 candidates are not dominated by any other candidate.
#' Rank 2 candidates are dominated only by Rank 1 candidates, etc.
#'
#' @param data A data.frame of metrics.
#' @param specs A named list where names match data columns, and values are
#'        1 (Maximize) or -1 (Minimize).
#' @return An integer vector of ranks.
#' @keywords internal
calculate_pareto_rank <- function(data, specs) {
  n <- nrow(data)
  if (n == 0) return(integer(0))

  # Normalize data directions so we always Maximize (multiply Minimize cols by -1)
  # matrix conversion for speed
  mat <- as.matrix(data[, names(specs)])
  for (col in names(specs)) {
    if (specs[[col]] == -1) {
      mat[, col] <- -mat[, col]
    }
  }

  ranks <- rep(NA, n)
  current_rank <- 1
  remaining_ids <- seq_len(n)

  while (length(remaining_ids) > 0) {
    # subset matrix for speed
    sub_mat <- mat[remaining_ids, , drop = FALSE]
    n_rem <- length(remaining_ids)
    is_dominated <- rep(FALSE, n_rem)

    # Compare every item i against every other item j
    # Optimization: In R, loops are slow, but for N ~ 20-100 (typical templates),
    # N^2 is trivial (10k ops).
    # A strictly dominated item is WORSE in at least one dimension
    # and NO BETTER in any dimension.

    for (i in 1:n_rem) {
      # We want to see if item i is dominated by ANY item j
      # If we find one j that dominates i, we stop checking i

      # Vectorized check: Compare row i against all rows
      # diffs > 0 means i is better
      # diffs < 0 means i is worse
      diffs <- t(t(sub_mat) - sub_mat[i, ])

      # i is dominated by j if:
      # j is never worse than i (min(diff) >= 0) AND
      # j is strictly better in at least one (max(diff) > 0)
      # Note: diffs is (j - i). So if j > i, diff > 0.

      # Check if ANY row j dominates i
      better_or_equal <- rowSums(diffs < 0) == 0 # i is never better than j (j >= i for all cols)
      strictly_better <- rowSums(diffs > 0) > 0  # j is strictly better in at least one

      if (any(better_or_equal & strictly_better)) {
        is_dominated[i] <- TRUE
      }
    }

    # Assign rank to non-dominated items
    frontier_indices <- remaining_ids[!is_dominated]
    ranks[frontier_indices] <- current_rank

    # Remove them for next iteration
    remaining_ids <- remaining_ids[is_dominated]
    current_rank <- current_rank + 1
  }

  return(ranks)
}

#' Export Design Results
#' @description Writes all final output files (FASTA, CSV, GFF3, GenBank).
#' @keywords internal
#'
export_design_results <- function(
    output_dir, design_name,
    variants_genomic,
    var_data, guides, repair_template,
    genome, txdb, edit_region,
    primers, all_probes, optimization_scheme) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # create large region of 1k BP around the edit window
  large_window <- promoters(
    edit_region, 1000, 1000 + width(edit_region),
    use.names = FALSE)
  large_window_seq <- getSeq(genome, large_window)[[1]]
  strand(large_window) <- "*"
  # Write FASTA file
  chrom <- DNAStringSet(large_window_seq)
  names(chrom) <- as.character(large_window)
  Biostrings::writeXStringSet(
    chrom, file.path(output_dir, paste0(design_name, ".fa")))

  if (length(guides) > 0) {
    guides$is_disabled_by_default <- NULL
    guides <- guides[, !apply(guides, 2, function(x) all(is.na(x))), drop = FALSE]
  }

  # Write files for each component using the helper function
  # all_probes is dt, remove seqnames, start, end, width, strand
  # var_data id long GRanges() (for each guide separately, repeated info)
  # primers is dt, nothing changed here, should work with current function
  write_component_files(output_dir, design_name, variants_genomic, "query")
  write_component_files(output_dir, design_name, guides, "guides")

  if (length(repair_template) > 0) {

    # Ensure columns exist (defaults for sorting safety)
    if (is.null(repair_template$any_overlaps_noncoding)) repair_template$any_overlaps_noncoding <- FALSE
    if (is.null(repair_template$pam_disrupted_count)) repair_template$pam_disrupted_count <- 0
    if (is.null(repair_template$seed_disrupted_count)) repair_template$seed_disrupted_count <- 0
    if (is.null(repair_template$total_disruption_count)) repair_template$total_disruption_count <- 0
    if (is.null(repair_template$total_cadd)) repair_template$total_cadd <- 999
    # Assumption: total_snp_quality_score is the sum/max of safety tiers (Lower is Better)
    if (is.null(repair_template$total_snp_quality_score)) repair_template$total_snp_quality_score <- 999
    if (is.null(repair_template$is_ag_risky)) repair_template$is_ag_risky <- FALSE
    repair_template$n_snvs <- lengths(strsplit(repair_template$snvs_introduced, ";"))

    # --- BINNING ---
    # New schema: 0 (best) to 5 (worst)
    # Base bin depends on PAM + total disruption.
    pam_disrupted <- repair_template$pam_disrupted_count >= 1
    total_disruption <- repair_template$total_disruption_count

    base_disruption_bin <- rep(5L, length(repair_template))
    has_useful_disruption <- total_disruption > 0

    base_disruption_bin[has_useful_disruption & pam_disrupted & total_disruption >= 2] <- 0L
    base_disruption_bin[has_useful_disruption & pam_disrupted & total_disruption == 1] <- 1L
    base_disruption_bin[has_useful_disruption & !pam_disrupted & total_disruption >= 3] <- 2L
    base_disruption_bin[has_useful_disruption & !pam_disrupted & total_disruption == 2] <- 3L
    base_disruption_bin[has_useful_disruption & !pam_disrupted & total_disruption == 1] <- 4L

    unsafe_penalty <- if (optimization_scheme == "disruption_first") {
      rep(0L, length(repair_template))
    } else if ("unsafe_snv_count" %in% names(mcols(repair_template))) {
      penalty <- repair_template$unsafe_snv_count
      penalty[is.na(penalty)] <- 0L
      as.integer(penalty)
    } else {
      rep(0L, length(repair_template))
    }
    repair_template$disruption_bin <- pmin(5L, base_disruption_bin + unsafe_penalty)

    # --- HIERARCHICAL SORTING ---
    ordering <- switch(
      optimization_scheme,

      # STRATEGY 1: BALANCED (Recommended)
      # Logic:
      # 1. Veto Tier 5 structural risks immediately.
      # 2. Prioritize Efficacy (Bin 1 > 2 > 3).
      # 3. Within the same efficacy bin, prioritize Safety (lower quality score).
      # 4. If efficacy and safety are equal, prioritize Parsimony (fewer SNVs).
      "balanced" = order(
        repair_template$any_overlaps_noncoding,    # Hard Safety Constraint
        repair_template$disruption_bin,            # Primary: Efficacy Confidence
        repair_template$total_snp_quality_score,   # Secondary: Safety
        repair_template$n_snvs,                    # Tertiary: Parsimony
        repair_template$is_ag_risky                # Tie-breaker
      ),

      # STRATEGY 2: DISRUPTION FIRST
      # Logic: Get the best Disruption Bin possible.
      # Use Parsimony as the secondary sort (don't care about subtle safety/CADD).
      "disruption_first" = order(
        repair_template$disruption_bin,
        repair_template$n_snvs,
        -repair_template$total_disruption_count     # More cuts is better here
      ),

      # STRATEGY 3: SAFETY FIRST
      # Logic: Safety is King. Only consider disruption if safety scores are identical.
      "safety_first" = order(
        repair_template$any_overlaps_noncoding,
        repair_template$total_snp_quality_score,
        repair_template$disruption_bin,
        repair_template$n_snvs
      ),

      # Fallback
      order(repair_template$disruption_bin, repair_template$n_snvs)
    )

    repair_template <- repair_template[ordering]

    # Reorder columns to make the CSV more readable for the user
    # Add disruption_bin to the output so users understand the sorting
    priority_cols <- c(
      "coords",
      "snvs_introduced",
      "n_snvs",
      "disruption_bin",
      "pam_disrupted_count",
      "seed_disrupted_count",
      "total_disruption_count",
      "aln_guide",
      "aln_template",
      "any_overlaps_noncoding",
      "total_cadd",
      "is_ag_risky",
      "total_snp_quality_score",
      "sequence"
    )

    # Identify columns that actually exist in the object
    existing_priority <- intersect(priority_cols, names(mcols(repair_template)))
    remaining_cols <- setdiff(names(mcols(repair_template)), existing_priority)

    # Reconstruct mcols in desired order
    mcols(repair_template) <- mcols(repair_template)[, c(existing_priority, remaining_cols)]
    # Internal helper, not part of exported templates schema
    repair_template$unsafe_snv_count <- NULL
    write_component_files(output_dir, design_name, as.data.frame(repair_template), "templates")
  }
  rownames(all_probes) <- all_probes$names
  write_component_files(output_dir, design_name, all_probes, "probes")

  if (!is.null(primers) && nrow(primers) > 0) { # Check nrow for data.frame
    # Calculate relative coordinates within the large_window window
    primers$LEFT_coords <- paste0(
      primers$chrom, ":", primers$LEFT_genomic_start, "-",
      as.numeric(primers$LEFT_genomic_start) + as.numeric(primers$LEFT_size))
    primers$RIGHT_coords <- paste0(
      primers$chrom, ":", primers$RIGHT_genomic_start, "-",
      as.numeric(primers$RIGHT_genomic_start) + as.numeric(primers$RIGHT_size))
    primers$RIGHT_genomic_start <- NULL
    primers$LEFT_genomic_start <- NULL
    primers$chrom <- NULL
    front <- intersect(c("names", "RIGHT_coords", "LEFT_coords"), colnames(primers))
    primers <- primers[, c(front, setdiff(names(primers), front)), drop = FALSE]

    write.table(
      as.data.frame(primers),
      file = file.path(output_dir, paste0(design_name, "_primers_1based.csv")),
      quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  }

  # We move the complex annotations to separate tables
  if (length(var_data) > 0) {
    # deduplicate
    var_names <- sapply(strsplit(names(var_data), " "), `[[`, 1)
    var_data <- var_data[!duplicated(var_names)]
    names(var_data) <- var_names[!duplicated(var_names)]
    var_data$guide_name <- NULL
    var_data$position_in_guide <- NULL
    var_data$ag_impact_score <- NULL
    var_data$ag_has_data <- NULL
    var_data$cadd_imputed <- NULL
    var_data$disruption_tier <- NULL
    var_data$priority_group <- NULL
    var_data$alphagenome_SPLICE_SITES <- NULL
    var_data$alphagenome_SPLICE_SITE_USAGE <- NULL
    var_data$alphagenome_SPLICE_JUNCTIONS <- NULL
    var_data_cds <- as.data.frame(var_data$CDS)
    if (nrow(var_data_cds) > 0) var_data_cds$group_name <- names(var_data)[var_data_cds$group]
    var_data_dbSNP <- as.data.frame(var_data$dbSNP)
    if (nrow(var_data_dbSNP) > 0) var_data_dbSNP$group_name <- names(var_data)[var_data_dbSNP$group]
    var_data_noncoding <- as.data.frame(var_data$noncoding)
    if (nrow(var_data_noncoding) > 0) var_data_noncoding$group_name <-
      names(var_data)[var_data_noncoding$group]
    var_data$CDS <- NULL
    var_data$dbSNP <- NULL
    var_data$noncoding <- NULL
    write_component_files(output_dir, design_name, var_data, "variants")

    if (!isEmpty(var_data_cds)) {
      write.table(var_data_cds, file = file.path(
        output_dir, paste0(design_name, "_cds_annotations_1based.csv")),
        quote = F, sep = ",", row.names = F, col.names = T)

      # Potentially unused code - TODO delete it?
      # unique_tx_ids <- unique(var_data_cds$tx_id)
      # if (length(unique_tx_ids) > 0) {
      #   cds <- suppressWarnings(cdsBy(txdb, by = "tx", use.names = T))
      #   cds <- cds[names(cds) %in% unique_tx_ids]
      #   cds <- suppressWarnings(restrict(cds, start = start(large_window), end = end(large_window), keep.all.ranges = TRUE))
      #   cds <- unlist(cds, use.names = T)
      #   cds <- cds[width(cds) > 0]
      #
      #   if (length(cds) > 0) {
      #     cds_on_large_window <- pmapToTranscripts(cds, large_window)
      #     write_genbank_custom(
      #       locus_name = as.character(large_window),
      #       gseq = large_window_seq,
      #       cds_on_tx = ranges(cds_on_large_window),
      #       aa_cds_seq = AAString(""), # Placeholder for AA sequence
      #       organism_name = organism(genome),
      #       output_file = file.path(output_dir, paste0(design_name, "_cds.gbk")))
      #   }
      #}
    }

    if (!isEmpty(var_data_dbSNP)) {
      write.table(var_data_dbSNP, file = file.path(
        output_dir, paste0(design_name, "_dbSNP_annotations_1based.csv")),
        quote = F, sep = ",", row.names = F, col.names = T)
    }

    if (!isEmpty(var_data_noncoding)) {
      write.table(var_data_noncoding, file = file.path(
        output_dir, paste0(design_name, "_noncoding_annotations_1based.csv")),
        quote = F, sep = ",", row.names = F, col.names = T)
    }
  }
}
