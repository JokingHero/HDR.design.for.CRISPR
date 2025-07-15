#' This new helper function abstracts the creation of a single template and its probes.
#' It will be used inside the switch statement to reduce code duplication.
#' @keywords internal
create_template_and_probes <- function(muts, # The list of mutations from a 'get_combinations' function
                                       guide_name,
                                       origin_mutation,
                                       template_range,
                                       genomic_seq,
                                       mutation_name,
                                       do_probes,
                                       probe_params) {

  mutated_seq <- replaceAt(genomic_seq,
                           at = ranges(muts[[1]]), value = DNAStringSet(muts[[1]]$replacement))

  rt <- GRanges(
    seqnames = mutation_name,
    ranges = template_range,
    strand = "+",
    replacement = as.character(extractAt(mutated_seq, template_range)),
    pam_disrupted = muts[[2]],
    guide_disrupted = muts[[3]],
    cadd = muts[[4]],
    overlaps_something = muts[[5]],
    snp_quality = -1 * muts[[6]]
  )
  # Create a unique name for the template based on the guide and mutations
  temp_name <- paste0("Template_", guide_name, "_Mut_", paste0(names(muts[[1]]), collapse = "_"))
  names(rt) <- temp_name

  hdr_probes <- GRanges()
  if (do_probes) {
    # design_probes takes a series of arguments, we use do.call to pass them as a list
    candidates <- do.call(design_probes, c(list(
      mutation_name = mutation_name,
      s = mutated_seq,
      origin_mut_start = start(origin_mutation)),
      probe_params))

    # Select the best probes to cover the introduced mutations
    hdr_probes <- select_probes(c(muts[[1]], origin_mutation), candidates, temp_name)
    names(hdr_probes) <- paste0("HDR probe ", seq_along(hdr_probes), " for ", temp_name)
  }

  list(template = rt, probes = hdr_probes)
}


#' Internal workhorse function for the HDR design pipeline.
#' @param params A list of all parameters passed from the user-facing function.
#' @keywords internal
design_hdr_internal <- function(params) {
  # Use with() to make parameters directly accessible without `params$` prefix
  with(params, {
    set.seed(seed) # Ensure reproducible randomness

    # --- 1. Initial Data Setup ---
    message("Setting up design environment...")
    env_data <- setup_design_environment(
      ensemble_transcript_id, mutation_loci, mutation_original, annotation, genome)

    message("Preparing candidate synonymous mutations...")
    mut_data <- prepare_candidate_snps(
      env_data, positions_to_mutate, mutation_loci, mutation_name,
      intron_bp, exon_bp, clinvar, snps, cadd, annotation, genome)

    # --- 2. Find and Filter Guides ---
    message("Finding and scoring guides...")
    guides <- get_guides_and_scores(mut_data$origin_mutation, mutation_name, guide_distance,
                                    env_data$genomic_seq, scores = score_efficiency)

    # Universal guide filter logic
    if (filter_to_guide != "") {
      found_guide <- guides$original == toupper(filter_to_guide)
      if (!any(found_guide)) {
        print(as.data.frame(guides))
        stop("Can't find the specified guide. See above for the guides we found.")
      }
      guides <- guides[found_guide]
      message(paste("Filtered to", length(guides), "guide(s)."))
    }
    pams <- resize(flank(guides, width = 3, start = FALSE), width = 2, fix = "end")

    # --- 3. Generate Templates based on Strategy ---
    message(paste("Generating templates with strategy:", strategy))
    repair_template <- GRanges()
    probes_ <- GRanges()
    template_range <- promoters(ranges(mut_data$origin_mutation),
                                upstream = template_upstream,
                                downstream = template_downstream)

    # Define common parameters for probe design
    probe_params <- list(
      st = start(mut_data$origin_mutation) + min(positions_to_mutate),
      sp = start(mut_data$origin_mutation) + max(positions_to_mutate),
      len_min = 19, len_max = 25, tmin = 50, tmax = 65)

    # The core switch based on the selected strategy
    switch(strategy,
           "optimal_per_guide" = {
             for (i in seq_along(guides)) {
               muts <- get_combinations_of_mutations_for_guide(
                 mut_data$mutations, maximum_mutations_per_template, pams[i], guides[i])

               if (is.null(muts) || length(muts[[1]]) == 0) {
                 message("Could not design synonymous SNPs for guide: ", names(guides[i]))
                 next
               }

               result <- create_template_and_probes(
                 muts, names(guides)[i], mut_data$origin_mutation,
                 template_range, env_data$genomic_seq[[1]],
                 mutation_name, do_probes, probe_params)
               repair_template <- c(repair_template, result$template)
               probes_ <- c(probes_, result$probes)
             }
           },
           "optimal_for_all" = {
             muts <- get_combinations_of_mutations_for_guides(
               mut_data$mutations, maximum_mutations_per_template, pams, guides)

             if (!is.null(muts)) {
               result <- create_template_and_probes(
                 muts, "All_Guides", mut_data$origin_mutation,
                 template_range, env_data$genomic_seq[[1]],
                 mutation_name, do_probes, probe_params)
               repair_template <- c(repair_template, result$template)
               probes_ <- c(probes_, result$probes)
             }
           },
           "all_per_guide" = {
             for (i in seq_along(guides)) {
               message("Working on guide: ", i, " (", guides$original[i], ")")
               for (mpt in 0:maximum_mutations_per_template) {

                 # Explicitly handle the case of zero synonymous SNPs to create a baseline template
                 # with only the original mutation corrected.
                 if (mpt == 0) {
                   # This creates a "muts" object with no mutations but correct structure
                   muts <- list(GRanges(), 0, 0, 0, FALSE, 0)
                   result <- create_template_and_probes(
                     muts, names(guides)[i], mut_data$origin_mutation,
                     template_range, env_data$genomic_seq[[1]],
                     mutation_name, do_probes, probe_params)
                   names(result$template) <- paste0("Template_", names(guides[i]), "_NoMut")
                   repair_template <- c(repair_template, result$template)
                   probes_ <- c(probes_, result$probes)
                   next # continue to next mpt
                 }

                 all_mut_combinations <- get_all_combinations_of_mutations_for_guide(
                   mut_data$mutations, mpt, pams[i], guides[i])

                 if (is.null(all_mut_combinations) || length(all_mut_combinations) == 0) {
                   message("Could not design synonymous SNPs for guide:
                           ", guides[i], " ", guides[i]$original," with ",
                           mpt, " mutations per template")
                   next
                 }

                 for (mut_combination in all_mut_combinations) {
                   result <- create_template_and_probes(
                     mut_combination, names(guides)[i], mut_data$origin_mutation,
                     template_range, env_data$genomic_seq[[1]],
                     mutation_name, do_probes, probe_params)
                   repair_template <- c(repair_template, result$template)
                   probes_ <- c(probes_, result$probes)
                 }
               }
             }
           },
           stop("Invalid 'strategy' specified. Must be one of 'optimal_per_guide', 'optimal_for_all', or 'all_per_guide'.")
    )

    # --- 4. Design Common Probes and Primers ---
    if (do_probes) {
      # Reference probes - outside of the SNP area
      # Length/Tm= 20-25 bp, 60+/-1 oC
      ref_seq <- IRanges(start = start(origin_mutation), width = 1) + 120
      ref_seq <- setdiff(
        ref_seq, IRanges(start = start(origin_mutation), width = 1) +
          max(abs(positions_to_mutate)))
      # two chunks that can be used for Ref probes
      rp <- design_probes(
        mutation_name, start(ref_seq[1]), end(ref_seq[1]), genomic_seq[[1]],
        origin_mut_start = start(origin_mutation))
      rp <- rp[which.max(rp$GC)]
      rp2 <- design_probes(
        mutation_name, start(ref_seq[2]), end(ref_seq[2]), genomic_seq[[1]],
        origin_mut_start = start(origin_mutation))
      rp2 <- rp2[which.max(rp2$GC)]
      ref_probes <- c(rp, rp2)
      names(ref_probes) <- paste0("Ref probe ", seq_along(ref_probes))

      # NHEJ probe - ~20 bp and have Tm values 58-60oC
      # contain the original mutated sequence &
      # no SNPs and no correction of the mutation
      nhej_probes <- design_probes(
        mutation_name,
        start(origin_mutation) + min(positions_to_mutate),
        start(origin_mutation) + max(positions_to_mutate),
        genomic_seq[[1]],
        len_min = 19, len_max = 25, tmin = 50, tmax = 65,
        origin_mut_start = start(origin_mutation))
      nhej_probes <- select_probes(origin_mutation, nhej_probes, "NHEJ origin mutation probe")
      names(nhej_probes) <- paste0("NHEJ probe ", seq_along(nhej_probes))

      probes_ <- c(probes_, ref_probes, nhej_probes)
    }

    primers <- if (primer3 != "") {
      design_primers(
        guide_distance, positions_to_mutate, env_data$genomic_seq,
        mut_data$origin_mutation, env_data$mut_genomic, genome, primer3)
    } else {
      NULL
    }

    # --- 5. Export All Results ---
    message("Exporting results...")
    export_design_results(
      output_dir, mutation_name, env_data$genomic_seq,
      mut_data$origin_mutation, mut_data$mutations, guides, repair_template,
      genome = genome,
      cds_on_tx = env_data$cds_on_tx,
      aa_cds_seq = env_data$aa_cds_seq,
      primers = primers, probes_ = probes_,
      score_efficiency = score_efficiency, one_for_all = (strategy == "optimal_for_all"))

    message("Design process complete.")
    return(invisible())
  })
}


#' @title Design HDR Repair Templates for CRISPR
#'
#' @description A comprehensive function to design HDR repair templates. It can generate templates
#' using different strategies: finding one optimal template for each guide, finding one optimal
#' template for all guides, or generating all possible template combinations.
#'
#' @param strategy The template generation strategy. Must be one of:
#' \itemize{
#'   \item \code{"optimal_per_guide"}: (Default) Designs one optimized template for each guide. Templates are optimized by selecting synonymous SNPs to disrupt guides while minimizing other functional impacts.
#'   \item \code{"optimal_for_all"}: Designs a single optimized template that works best across all identified guides.
#'   \item \code{"all_per_guide"}: Generates all possible template combinations for each guide, for a range of mutation counts from 0 to \code{maximum_mutations_per_template}.
#' }
#' @param maximum_mutations_per_template The maximum number of synonymous SNPs to introduce per repair template.
#' @param filter_to_guide Optional. A 20bp character string of a specific guide sequence. If provided, the design process will be restricted to only this guide.
#' @param ensemble_transcript_id This has to be an Ensembl transcript ID, e.g., "ENST00000307851.9".
#' @param mutation_loci Position on the CDS that is mutated, e.g., 245.
#' @param mutation_original The original base on the genome, e.g., "A".
#' @param mutation_replacement The desired mutated base, e.g., "G".
#' @param mutation_name A short name for this mutation, e.g., "tyr82cys".
#' @param output_dir Path to the directory where output files will be saved.
#' @param annotation File path to the genome annotation file (GFF3 format).
#' @param genome A \code{BSgenome} object for your genome, e.g., \code{BSgenome.Hsapiens.UCSC.hg38}.
#' @param guide_distance Window around the mutation to search for guides (default: 27).
#' @param positions_to_mutate A vector of positions relative to the mutation to consider for introducing synonymous SNPs (default: -19:19).
#' @param template_upstream How many bp upstream of the mutation to include in the repair template (default: 59).
#' @param template_downstream How many bp downstream of the mutation to include in the repair template (default: 61).
#' @param intron_bp Number of intronic bases near splice sites to exclude from synonymous SNPs (default: 6).
#' @param exon_bp Number of exonic bases near splice sites to exclude from synonymous SNPs (default: 3).
#' @param snps Optional. An object like \code{SNPlocs.Hsapiens.dbSNP155.GRCh38}.
#' @param clinvar Optional. File path to a ClinVar VCF file.
#' @param cadd Optional. A function to retrieve CADD scores, e.g., `getGScores("cadd.v1.6.hg38")`.
#' @param score_efficiency Logical. If `TRUE`, score guides using models from the `crisprScore` package.
#' @param seed An integer for setting the random seed to ensure reproducibility.
#' @param do_probes Logical. If `TRUE`, design qPCR probes for HDR, NHEJ, and Reference.
#' @param primer3 Optional. Full path to the `primer3_core` executable to design PCR primers.
#'
#' @return This function does not return a value. It writes all output files to the specified `output_dir`.
#' @import Biostrings GenomicFeatures GenomicRanges SummarizedExperiment IRanges BSgenome BSgenome.Hsapiens.UCSC.hg38 VariantAnnotation GenomeInfoDb
#' @importFrom utils write.table combn
#' @export
design_hdr <- function(
    ensemble_transcript_id,
    mutation_loci,
    mutation_original,
    mutation_replacement,
    mutation_name,
    output_dir,
    annotation,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    strategy = c("optimal_per_guide", "optimal_for_all", "all_per_guide"),
    maximum_mutations_per_template = 3,
    filter_to_guide = "",
    guide_distance = 10 + 17,
    positions_to_mutate = -19:19,
    template_upstream = 59,
    template_downstream = 61,
    seed = 42,
    intron_bp = 6,
    exon_bp = 3,
    snps = NULL,
    clinvar = NULL,
    cadd = NULL,
    score_efficiency = FALSE,
    do_probes = TRUE,
    primer3 = ""
) {
  # Match the strategy argument to ensure it's one of the allowed options
  strategy <- match.arg(strategy)

  # Capture all function arguments into a list to pass to the internal engine
  params <- as.list(environment())

  # Call the internal workhorse function
  design_hdr_internal(params)
}
