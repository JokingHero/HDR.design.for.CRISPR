#' @title Design HDR Repair Templates for CRISPR
#'
#' @description A comprehensive function to design HDR repair templates. It can generate templates
#' using different strategies: finding one optimal template for each guide, finding one optimal
#' template for all guides, or generating all possible template combinations.
#'
#' @param design_name A short name for this set of mutations, e.g., "tyr82cys".0,
#' @param chrom Chromosome of the variant to correct. e.g. chr1
#' @param variant_start Position on the chromosome for the variants (its on + strand)
#' @param variant_end Position ends of the variants on the chromosome
#' @param REF The original bases on the genome, e.g., "A".
#' @param ALT The desired mutated bases, e.g., "G".
#' @param ALT_on_genome A vector of T/F for each variant. Guides should match your target cell sequences. Select false when targeting wild type.
#' @param ALT_on_templates A vector of T/F for each variant. This should match your desired outcome. Select true when correcting alt variant.
#' @param output_dir Path to the directory where output files will be saved.
#' @param annotation File path to the genome annotation file (GFF3/GTF) or (.db/.sqlite) that can be loaded with `AnnotationDbi::loadDb`.
#' @param genome A \code{BSgenome} object for your genome, e.g., \code{BSgenome.Hsapiens.UCSC.hg38}.
#' @param optimization_scheme The optimization scheme for selecting synonymous SNPs. Must be one of:
#' \itemize{
#'   \item \code{"balanced"}: (Default) Balances safety, PAM disruption, and SNP quality.
#'   \item \code{"safety_first"}: Prioritizes avoiding non-coding overlaps and using known benign SNPs.
#'   \item \code{"disruption_first"}: Prioritizes PAM disruption and guide disruption.
#' }
#' @param maximum_mutations_per_template The maximum number of synonymous SNPs to introduce per repair template.
#' @param filter_to_guide Optional. A 20bp character string of a specific guide sequence. If provided, the design process will be restricted to only this guide.
#' @param cut_distance_max Window around the variants to search for cut of the guides (default: 30).
#' @param template_size Size of the repair template (default: 120).
#' @param template_hr_arm_size Size of the repair template (default: 30).
#' @param seed An integer for setting the random seed to ensure reproducibility.
#' @param intron_bp Number of intronic bases near splice sites to exclude from synonymous SNPs (default: 6).
#' @param exon_bp Number of exonic bases near splice sites to exclude from synonymous SNPs (default: 3).
#' @param snps Optional. An object like \code{SNPlocs.Hsapiens.dbSNP155.GRCh38}.
#' @param clinvar Optional. File path to a ClinVar VCF file.
#' @param cadd Optional. A function to retrieve CADD scores, e.g., `getGScores("cadd.v1.6.hg38")`.
#' @param score_efficiency Logical. If `TRUE`, score guides using models from the `crisprScore` package.
#' @param do_probes Logical. If `TRUE`, design qPCR probes for HDR, NHEJ, and Reference.
#' @param primer3 Optional. Full path to the `primer3_core` executable to design PCR primers.
#' @param alphagenome_key Your key to the alphagenome service.
#' @param alphagenome_context Celltypes for which to filter alphagenome scores.
#' @param python_exec path to the python3 that has alphagenome installed.
#' @return This function does not return a value. It writes all output files to the specified `output_dir`.
#' @import Biostrings GenomicFeatures GenomicRanges SummarizedExperiment IRanges BSgenome BSgenome.Hsapiens.UCSC.hg38 VariantAnnotation GenomeInfoDb
#' @importFrom utils write.table combn
#' @export
#'
design_hdr <- function(
  design_name,
  chrom,
  variant_start,
  variant_end,
  REF,
  ALT,
  ALT_on_genome,
  ALT_on_templates,
  output_dir,
  annotation,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  optimization_scheme = "balanced",
  maximum_mutations_per_template = 3,
  filter_to_guide = "",
  cut_distance_max = 30,
  template_size = 120,
  template_hr_arm_size = 30,
  seed = 42,
  intron_bp = 6,
  exon_bp = 3,
  snps = NULL,
  clinvar = NULL,
  cadd = NULL,
  score_efficiency = FALSE,
  do_probes = TRUE,
  primer3 = "",
  alphagenome_key = "",
  alphagenome_context = "",
  python_exec = "python3"
) {
  set.seed(seed) # Ensure reproducible randomness

  stopifnot(length(variant_start) == length(variant_end))
  stopifnot(length(ALT_on_genome) == length(variant_start))
  stopifnot(length(ALT_on_templates) == length(variant_start))
  stopifnot(length(ALT) == length(variant_start))
  stopifnot(length(REF) == length(variant_start))

  optimization_scheme <- match.arg(optimization_scheme, c("balanced", "safety_first", "disruption_first"))
  set.seed(seed)

  # --- 1. Initial Data Setup ---
  variants_genomic <- GRanges(
    seqnames = chrom,
    ranges = IRanges(
      start = variant_start,
      end = variant_end
    ),
    strand = "+",
    REF = REF, ALT = ALT
  )
  variants_genomic <- normalize_variants(variants_genomic, genome)
  effective_ts <- template_size - 2 * template_hr_arm_size
  if (width(range(variants_genomic)) > effective_ts) {
    stop("Your variants are spaced too wide for this template size.")
  }

  txdb <- if (tools::file_ext(annotation) %in% c("gff3", "gtf")) {
    suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation))
  } else {
    suppressWarnings(AnnotationDbi::loadDb(annotation))
  }

  message("Finding and scoring guides...")
  guides <- get_guides_and_scores(
    variants_genomic, design_name,
    cut_distance_max, genome, ALT_on_genome, score_efficiency)

  # Universal guide filter logic
  if (filter_to_guide != "") {
    found_guide <- guides$original == toupper(filter_to_guide)
    if (!any(found_guide)) {
      stop(paste(
        "Can't find the specified guide. We found:\n",
        paste0(guides$original, collapse = "\n")
      ))
    }
    guides <- guides[found_guide]
  }

  message("Preparing candidate synonymous mutations...")
  # Calculate the net length change from variants destined for the template
  net_length_change <- if (any(ALT_on_templates)) {
    sum(nchar(variants_genomic$ALT[ALT_on_templates])) -
      sum(nchar(variants_genomic$REF[ALT_on_templates]))
  } else {
    0
  }

  # Define the genomic footprint size needed to achieve the final template_size
  genomic_footprint_size <- template_size - net_length_change
  if (genomic_footprint_size <= 0) {
    stop("The desired template_size is too small to accommodate the large insertions from your variants.")
  }
  variants_center_gr <- if (any(ALT_on_templates)) {
    range(variants_genomic[ALT_on_templates])
  } else {
    range(variants_genomic)
  }
  # This `edit_region` is our definitive genomic range for the final HDR template
  # after we inject variants it will have proper HDR template size
  edit_region <- resize(variants_center_gr, width = genomic_footprint_size, fix = "center")

  # Create the base HDR template sequence by extracting the genomic footprint and injecting variants
  source_genomic_seq <- suppressWarnings(getSeq(genome, edit_region))[[1]]
  variants_in_editw <- pmapToTranscripts(variants_genomic[ALT_on_templates], edit_region)
  mcols(variants_in_editw) <- mcols(variants_genomic[ALT_on_templates])
  names(variants_in_editw) <- names(variants_genomic[ALT_on_templates])

  preHDR_seq <- if (length(variants_in_editw) > 0) {
    replaceAt(source_genomic_seq,
              at = ranges(variants_in_editw),
              value = DNAStringSet(variants_in_editw$ALT)
    )
  } else {
    source_genomic_seq
  }

  # Align guides to the new template sequence to find their locations
  aln <- pwalign::pairwiseAlignment(
    DNAStringSet(guides$with_pam), preHDR_seq,
    type = "global-local",
    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
      match = 1, mismatch = 0, baseOnly = F
    ),
    gapOpening = 0,
    gapExtension = -1
  )

  # --- Create a detailed map of potential SNP locations linked to specific guides ---
  guides_on_template <- GRanges("target_seq", ranges = ranges(subject(aln)))
  names(guides_on_template) <- rownames(guides)
  tiled_guides_on_template_list <- tile(guides_on_template, width = 1)
  position_in_footprint <- unlist(lapply(
    tiled_guides_on_template_list,
    function(gr) seq_along(gr)
  ))
  tiled_guides_on_template <- unlist(tiled_guides_on_template_list,
                                     use.names = FALSE
  )
  is_plus_strand <- guides$strand[
    match(names(tiled_guides_on_template), rownames(guides))
  ] == "+"

  # --- Protect homology arms from synonymous mutations ---
  active_region_on_template <- IRanges(
    start = template_hr_arm_size + 1,
    end = nchar(preHDR_seq) - template_hr_arm_size
  )
  is_in_active_region <- overlapsAny(
    ranges(tiled_guides_on_template),
    active_region_on_template,
    type = "within"
  )
  tiled_guides_on_template <- tiled_guides_on_template[is_in_active_region]
  position_in_footprint <- position_in_footprint[is_in_active_region]
  is_plus_strand <- is_plus_strand[is_in_active_region]

  # TODO this is not necesairly stop, this might be that there are some
  # disabled guides! that are valid and we can output these
  if (length(tiled_guides_on_template) == 0) {
    stop("No potential SNP positions available within the active template region after excluding homology arms.")
  }

  position_in_guide <- integer(length(tiled_guides_on_template))
  position_in_guide[is_plus_strand] <- position_in_footprint[is_plus_strand]
  position_in_guide[!is_plus_strand] <- 20 - (position_in_footprint[!is_plus_strand] - 4)

  mcols(tiled_guides_on_template) <- NULL
  tiled_guides_on_template$guide_name <- names(tiled_guides_on_template)
  tiled_guides_on_template$position_in_guide <- position_in_guide

  # Map candidate SNP locations from template coordinates to genomic coordinates
  template_coord_map <- if (any(ALT_on_templates)) {
    build_variant_layout(variants_in_editw, nchar(source_genomic_seq))
  } else {
    GRanges("target_seq", IRanges(1, nchar(source_genomic_seq)),
            source = "genomic", origin_id = "genomic",
            origin_start = 1, origin_end = nchar(source_genomic_seq)
    )
  }

  names(tiled_guides_on_template) <- NULL
  candidate_snp_map <- remap_target_to_genomic(
    target = tiled_guides_on_template,
    coordinate_map = template_coord_map,
    window_genomic = edit_region,
    variants_genomic = variants_genomic[ALT_on_templates]
  )
  # This filters SNVs positions that are over potential Variants
  # because that has to be preserved!
  candidate_snp_map <- candidate_snp_map[
    !startsWith(candidate_snp_map$coords, "Variant"),
  ]

  variants_genomic_on_ts <- variants_genomic[ALT_on_templates]
  var_data <- prepare_candidate_snps(
    candidate_snp_map,
    annotation, txdb, genome,
    variants_genomic_on_ts,
    intron_bp, exon_bp, clinvar, snps, cadd,
    alphagenome_key, python_exec, alphagenome_context
  )

  # TODO: we could still have valid self deactivated guides here!!!
  # same as above
  if (length(var_data) == 0) {
    stop("No valid candidate SNPs could be generated for the active guides.")
  }
  var_data <- augment_var_data_with_scores(var_data, optimization_scheme)

  # --- Generate Templates based on Optimization Scheme ---
  message(paste("Generating templates with scheme:", optimization_scheme))
  repair_template <- GRanges()
  probes_ <- list()
  probe_params <- list(len_min = 19, len_max = 25, tmin = 50, tmax = 65)

  for (guide_id in rownames(guides)) {
    message("Working on guide: ", guide_id)
    current_guide <- guides[guide_id, ]
    guide_snps <- if (length(var_data) > 0) var_data[var_data$guide_name == guide_id] else GRanges()
    for (mpt in 0:maximum_mutations_per_template) {
      selected_muts <- GRanges()
      if (mpt > 0) {
        # Select optimal subset of size 'mpt'
        if (length(guide_snps) > 0) {
          selected_muts <- find_best_snps_for_guide(
            guide_snps, mpt, optimization_scheme
          )
          # we can drop extra (NGG_1) in names as its for this specific guide anyway
          names(selected_muts) <- sapply(strsplit(names(selected_muts), " "), `[[`, 1)
        }
        # If we couldn't find 'mpt' valid non-overlapping mutations, skip this iteration
        if (length(selected_muts) != mpt) next
      }

      # Generate Template and Calculate Scores
      design_id <- paste0(guide_id, "_with_", mpt, "_SNVs")
      result <- create_template_and_probes(
        selected_muts, variants_in_editw, variants_genomic_on_ts,
        design_id, edit_region, source_genomic_seq,
        current_guide,
        do_probes, probe_params
      )

      repair_template <- c(repair_template, result$template)
      probes_[[design_id]] <- result$probes
    }
  }

  # --- 4. Design Control Probes and Amplification Primers ---

  # Consolidate all probes into one GRanges object
  all_probes <- do.call(rbind, probes_)
  rownames(all_probes) <- NULL

  if (do_probes) {
    message("Designing NHEJ control probes...")
    # These probes are the probes that need to verify whether HDR is positively
    # integrated in negative fashion
    # they are on genomic sequence "ALT_on_genome"
    variants_for_guides_genomic <- variants_genomic[ALT_on_genome]
    variants_for_guides_in_editw <- pmapToTranscripts(variants_for_guides_genomic, edit_region)
    mcols(variants_for_guides_in_editw) <- mcols(variants_for_guides_genomic)
    names(variants_for_guides_in_editw) <- names(variants_for_guides_genomic)

    nhej_target_seq <- if (length(variants_for_guides_in_editw) > 0) {
      replaceAt(source_genomic_seq,
                at = ranges(variants_for_guides_in_editw),
                value = DNAStringSet(variants_for_guides_in_editw$ALT)
      )
    } else {
      source_genomic_seq
    }

    # We need a coordinate map for the nhej_target_seq
    nhej_coord_map <- build_variant_layout(
      variants_for_guides_in_editw, nchar(source_genomic_seq)
    )

    nhej_probe_candidates <- design_probes(
      s = nhej_target_seq,
      genomic_context = edit_region,
      coordinate_map = nhej_coord_map,
      variants_genomic = variants_for_guides_genomic,
      tmin = probe_params$tmin, tmax = probe_params$tmax,
      len_min = probe_params$len_min, len_max = probe_params$len_max
    )

    # 3. Select probes that specifically cover the loci
    # targeted for change by the HDR template.
    if (nrow(nhej_probe_candidates) > 0) {
      # but they should check all positions of ALL variants
      nhej_probes <- select_probes(
        muts_to_cover = variants_genomic,
        candidates = nhej_probe_candidates,
        temp_name = "NHEJ_control"
      )
      if (nrow(nhej_probes) > 0) {
        nhej_probes$names <- paste0("NHEJ_probe_", seq_len(nrow(nhej_probes)))
        rownames(nhej_probes) <- NULL
        all_probes <- rbind(all_probes, nhej_probes)
      }
    } else {
      warning("Could not generate any candidate probes for NHEJ detection.")
    }
  }

  primers <- if (primer3 != "") {
    message("Designing amplification primers...")
    # Define a larger region to search for primers, flanking the edit_region.
    primer_search_region <- promoters(edit_region, 150, 150 + width(edit_region))

    # Ensure the product must contain the full span of template variants.
    design_primers(
      template_range_extended = primer_search_region,
      template_range = edit_region,
      variant_with_allowed = variants_center_gr,
      genome = genome,
      primer3_path = primer3
    )
  } else {
    NULL
  }

  message("Exporting results...")
  variants_genomic$ALT_on_genome <- ALT_on_genome
  variants_genomic$ALT_on_templates <- ALT_on_templates
  export_design_results(
    output_dir, design_name,
    variants_genomic, var_data, guides, repair_template,
    genome, txdb, edit_region,
    primers = primers,
    all_probes = all_probes,
    optimization_scheme = optimization_scheme
  )

  message("Design process complete.")
  return(invisible())
}
