#' This new helper function abstracts the creation of a single template and its probes.
#' It will be used inside the switch statement to reduce code duplication.
#' @keywords internal
#'
create_template_and_probes <- function(muts,
                                       guide_name,
                                       template_range,
                                       template_ref,
                                       design_name,
                                       do_probes,
                                       probe_params){
  muts_on_temp <- pmapToTranscripts(muts$mutations, template_range)
  # Its very important to have IRanges for replacement in replaceAt
  # 0 width ranges or numeric vector is an insertion instead
  mutated_seq <- replaceAt(
    template_ref, at = IRanges(start(muts_on_temp), width = 1),
    value = DNAStringSet(muts$mutations$ALT))
  template_range$REF <- NULL
  template_range$ALT <- as.character(mutated_seq)
  template_range$pam_disrupted <- muts$pam_disrupted_count
  template_range$guide_disrupted <- muts$guide_disrupted_count
  template_range$cadd <- muts$total_cadd
  template_range$overlaps_noncoding <- muts$any_overlaps_noncoding
  template_range$snp_quality <- -1 * muts$total_compatibility_score

  # Create a unique name for the template based on the guide and mutations
  temp_name <- paste0("Template_", guide_name, "_Mut_",
                      paste0(names(muts$mutations), collapse = "; "))
  names(template_range) <- temp_name

  hdr_probes <- GRanges()
  if (do_probes) {
    # design_probes takes a series of arguments, we use do.call to pass them as a list
    candidates <- do.call(design_probes, c(list(
      s = mutated_seq,
      template_range = template_range),
      probe_params))

    # Select the best probes to cover the introduced mutations
    hdr_probes <- select_probes(muts$mutations, candidates, temp_name)
    names(hdr_probes) <- paste0("HDR probe ", seq_along(hdr_probes), " for ", temp_name)
  }

  list(template = template_range, probes = hdr_probes)
}

#' @title Design HDR Repair Templates for CRISPR
#'
#' @description A comprehensive function to design HDR repair templates. It can generate templates
#' using different strategies: finding one optimal template for each guide, finding one optimal
#' template for all guides, or generating all possible template combinations.
#'
#' @param design_name A short name for this mutation, e.g., "tyr82cys".0,
#' @param chrom Chromosome of the variant to correct. e.g. chr1
#' @param variant_start Position on the chromosome for the variant (its on + strand)
#' @param variant_end Position end of the variant on the chromosome
#' @param REF The original base on the genome, e.g., "A".
#' @param ALT The desired mutated base, e.g., "G".
#' @param ALT_on_guides Guides should match your target cell sequences. Select false when targetting wild type.
#' @param ALT_on_templates This should match your desired outcome. Select true when correcting alt variant.
#' @param output_dir Path to the directory where output files will be saved.
#' @param annotation File path to the genome annotation file (GFF3/GTF) or (.db/.sqlite) that can be loaded with `AnnotationDbi::loadDb`.
#' @param genome A \code{BSgenome} object for your genome, e.g., \code{BSgenome.Hsapiens.UCSC.hg38}.
#' @param strategy The template generation strategy. Must be one of:
#' \itemize{
#'   \item \code{"optimal_per_guide"}: (Default) Designs one optimized template for each guide. Templates are optimized by selecting synonymous SNPs to disrupt guides while minimizing other functional impacts.
#'   \item \code{"optimal_for_all"}: Designs a single optimized template that works best across all identified guides.
#'   \item \code{"all_per_guide"}: Generates all possible template combinations for each guide, for a range of mutation counts from 0 to \code{maximum_mutations_per_template}.
#' }
#' @param maximum_mutations_per_template The maximum number of synonymous SNPs to introduce per repair template.
#' @param filter_to_guide Optional. A 20bp character string of a specific guide sequence. If provided, the design process will be restricted to only this guide.
#' @param guide_distance Window around the mutation to search for guides (default: 27).
#' @param allowed_positions_upstream A vector of positions upstream to the variant to consider for introducing synonymous SNPs (default: 1:10).
#' @param allowed_positions_downstream A vector of positions downstream to the variant to consider for introducing synonymous SNPs (default: 1:10).
#' @param template_upstream How many bp upstream of the mutation to include in the repair template (default: 59).
#' @param template_downstream How many bp downstream of the mutation to include in the repair template (default: 61).
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
#' @param python_exec path to the python3 that has alphagenome installed.
#' @return This function does not return a value. It writes all output files to the specified `output_dir`.
#' @import Biostrings GenomicFeatures GenomicRanges SummarizedExperiment IRanges BSgenome BSgenome.Hsapiens.UCSC.hg38 VariantAnnotation GenomeInfoDb
#' @importFrom utils write.table combn
#' @export
design_hdr <- function(
    design_name,
    chrom,
    variant_start,
    variant_end,
    REF,
    ALT,
    ALT_on_guides,
    ALT_on_templates,
    output_dir,
    annotation,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    strategy,
    maximum_mutations_per_template = 3,
    filter_to_guide = "",
    guide_distance = 10 + 17,
    allowed_positions_upstream = 1:10,
    allowed_positions_downstream = 1:10,
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
    primer3 = "",
    alphagenome_key = "",
    python_exec = "python3"
) {
  strategy <- match.arg(strategy, c("optimal_per_guide", "optimal_for_all", "all_per_guide"))
  set.seed(seed)
  variant_genomic <- GRanges(seqnames = chrom,
                             ranges = IRanges(start = variant_start,
                                              end = variant_end),
                             strand = "+",
                             REF = REF, ALT = ALT)
  txdb <- if (tools::file_ext(annotation) %in% c("gff3", "gtf")) {
    suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation))
  } else {
    suppressWarnings(AnnotationDbi::loadDb(annotation))
  }

  message("Preparing candidate synonymous mutations...")
  var_data <- prepare_candidate_snps(
    annotation, txdb, genome,
    variant_genomic, allowed_positions_upstream, allowed_positions_downstream,
    intron_bp, exon_bp, clinvar, snps, cadd,
    alphagenome_key, python_exec)

  message("Finding and scoring guides...")
  guides <- get_guides_and_scores_refactored(
    variant_genomic, design_name, guide_distance, genome, ALT_on_guides, score_efficiency)
  if (filter_to_guide != "") {
    found_guide <- guides$original == toupper(filter_to_guide)
    if (!any(found_guide)) {
      print(as.data.frame(guides))
      stop("Can't find the specified guide. See above for the guides we found.")
    }
    guides <- guides[found_guide]
  }
  pams <- resize(flank(guides, width = 3, start = FALSE), width = 2, fix = "end")

  message(paste("Generating templates with strategy:", strategy))
  repair_template <- GRanges()
  probes_ <- GRanges()
  template_range <- promoters(variant_genomic,
                              upstream = template_upstream,
                              downstream = template_downstream)
  strand(template_range) <- "+"
  template_ref <- getSeq(genome, template_range)[[1]]
  variant_temp <- pmapToTranscripts(variant_genomic, template_range)
  template_alt <- replaceAt(
    template_ref, at = ranges(variant_temp),
    value = DNAStringSet(variant_genomic$ALT))
  variant_with_allowed <- promoters(
    variant_genomic,
    upstream = max(allowed_positions_upstream),
    downstream = max(allowed_positions_downstream))

  probe_params <- list(len_min = 19, len_max = 25, tmin = 50, tmax = 65)
  switch(strategy,
         "optimal_per_guide" = {
           for (i in seq_along(guides)) {
             muts <- find_mutation_combinations(
               var_data,  guides[i], pams[i], maximum_mutations_per_template, "best_single")

             if (is.null(muts) || length(muts[[1]]) == 0) {
               message("Could not design synonymous SNPs for guide: ", names(guides[i]))
               next
             }

             result <- create_template_and_probes(
               muts, names(guides)[i],
               template_range, if (ALT_on_templates) template_alt else template_ref,
               design_name, do_probes, probe_params)
             repair_template <- c(repair_template, result$template)
             probes_ <- c(probes_, result$probes)
           }
         },
         "optimal_for_all" = {
           muts <- find_mutation_combinations(
             var_data, guides, pams, maximum_mutations_per_template, "best_global")

           if (!is.null(muts)) {
             result <- create_template_and_probes(
               muts, "All_Guides",
               template_range, if (ALT_on_templates) template_alt else template_ref,
               design_name, do_probes, probe_params)
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
                 muts <- list(
                   mutations = variant_genomic,
                   pam_disrupted_count = 0,
                   guide_disrupted_count = 0,
                   total_cadd = 0,
                   any_overlaps_noncoding = FALSE,
                   total_compatibility_score = 0)
                 muts$mutations$ALT <- muts$mutations$REF
                 # this is "corrected" variant, which is the same as
                 # on template_ref, but we will find probe over this variant
                 result <- create_template_and_probes(
                   muts, names(guides)[i],
                   template_range, if (ALT_on_templates) template_alt else template_ref,
                   design_name, do_probes, probe_params)
                 names(result$template) <- paste0("Template_", names(guides[i]), "_NoMut")
                 repair_template <- c(repair_template, result$template)
                 probes_ <- c(probes_, result$probes)
                 next # continue to next mpt
               }

               all_mut_combinations <- find_mutation_combinations(
                 var_data, guides[i], pams[i], mpt, "all_valid")

               if (is.null(all_mut_combinations) || length(all_mut_combinations) == 0) {
                 message("Could not design synonymous SNPs for guide:
                           ", guides[i], " ", guides[i]$original," with ",
                         mpt, " mutations per template")
                 next
               }

               for (mut_combination in all_mut_combinations) {
                 result <- create_template_and_probes(
                   mut_combination, names(guides)[i],
                   template_range, if (ALT_on_templates) template_alt else template_ref,
                   design_name, do_probes, probe_params)
                 repair_template <- c(repair_template, result$template)
                 probes_ <- c(probes_, result$probes)
               }
             }
           }
         },
         stop("Invalid 'strategy' specified. Must be one of 'optimal_per_guide', 'optimal_for_all', or 'all_per_guide'.")
  )

  if (do_probes) {
    # Reference probes - outside of the SNP area
    # Length/Tm= 20-25 bp, 60+/-1 oC
    upstream_template_range <- template_range
    end(upstream_template_range) <- GenomicRanges::start(variant_with_allowed)
    start(upstream_template_range) <- start(upstream_template_range) - 50
    downstream_template_range <- template_range
    start(downstream_template_range) <- GenomicRanges::end(variant_with_allowed)
    end(downstream_template_range) <- end(downstream_template_range) + 50

    # NHEJ probe - ~20 bp and have Tm values 58-60oC
    # contain the original mutated sequence &
    # no SNPs and no correction of the mutation
    nhej_probes <- design_probes(if (ALT_on_templates) template_ref else template_alt, template_range,
                                 len_min = 19, len_max = 25, tmin = 50, tmax = 65)
    nhej_probes <- select_probes(variant_genomic, nhej_probes, "NHEJ origin mutation probe")
    names(nhej_probes) <- paste0("NHEJ probe ", seq_along(nhej_probes))

    probes_ <- c(probes_, nhej_probes)
  }

  primers <- if (primer3 != "") {
    # we need to make now target_seq that is bigger than HDR template
    # we also want at least one primer on the genome outside of HDR template
    template_range_extended <- promoters(template_range, 150, 150 + width(template_range))
    design_primers(
      template_range_extended, template_range, variant_with_allowed, genome, primer3)
  } else {
    GRanges()
  }

  message("Exporting results...")
  export_design_results(
    output_dir, design_name,
    variant_genomic, var_data, guides, repair_template,
    genome, txdb,
    primers = primers, probes_ = probes_,
    score_efficiency = score_efficiency,
    one_for_all = (strategy == "optimal_for_all"))

  message("Design process complete.")
  return(invisible())
}
