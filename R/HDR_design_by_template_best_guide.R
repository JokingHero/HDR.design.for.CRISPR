# rm(list = ls(all.names = TRUE))
# gc(reset = T)
#
# library(Biostrings)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(IRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(crisprScore)
#
# ensemble_transcript_id = "ENST00000399837.8"
# mutation_loci = 506
# mutation_original = "G"
# mutation_replacement = "A"
# mutation_name = "R169Q"
# output_dir = "~/"
# annotation = "/mnt/corsair/Projects/uib/CRIPSR/Oslo/Oline_seqeuncing_HDR_14_07_2023/gencode.v42.annotation.gff3.gz"
# genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# guide_distance = 10 + 17
# extension = 400
# positions_to_mutate = -19:19
# mutations_per_template = 4
# seed = 42

# ensemble_transcript_id <- "ENST00000374202.7"
# mutation_loci <- 172
# mutation_original <- "C"
# mutation_replacement <- "T"
# mutation_name <- "pro58sel"
#
# # HAV2RC with the mutation c.245A>G (p.Tyr82Cys).
# ensemble_transcript_id <- "ENST00000307851.9"
# mutation_loci <- 245
# mutation_original <- "A"
# mutation_replacement <- "G"
# mutation_name <- "tyr82cys"

#' @title Design best CRISPR guide for fixing mutation
#'
#' @description Design CRISPR guide and HDR repair template so that it auto-inhibits the guide.
#'
#' @param ensemble_transcript_id This has to be ensemble transcript id e.g. ENST00000307851.9
#' @param mutation_loci Position on the gene that is mutated e.g. 245
#' @param mutation_original What was the original base on the genome? e.g. A
#' @param mutation_replacement What is the mutated base? e.g. G
#' @param mutation_name What is your name for this mutation e.g. tyr82cys
#' @param output_dir Where the files should be generated? Make sure you have writing permissions there.
#' @param annotation File path to the annotation file, a gff3.
#' @param genome BSgenome of your genome, compatible with annotation file, by default it is hg38.
#' @param guide_distance Window around which we should search for guides, relatively to the mutation loci. Defualt is 10 + 17bp.
#' @param extension How many bases upstream/downstream from the mutation loci should we include. Default is 400.
#' @param positions_to_mutate Which positions from the mutation are available for mutation. By default its -30:30
#' leaving 20bp on each side of the template for the in case of incomplete integration of the template. Also codon
#' occupied by the mutation is not included in the change.
#' @param mutations_per_template How many codons should be mutated per template.
#' @param seed Ensures reproducibility of the random parts of the pipeline.
#' @return writes files to the specified directory, might overwrite
#' @import Biostrings GenomicFeatures GenomicRanges IRanges BSgenome.Hsapiens.UCSC.hg38
#' @importFrom utils write.table
#' @export
#'
design_template_for_best_guide <-
  function(ensemble_transcript_id,
           mutation_loci,
           mutation_original,
           mutation_replacement,
           mutation_name,
           output_dir,
           annotation,
           genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
           guide_distance = 10 + 17,
           extension = 400,
           positions_to_mutate = -19:19,
           mutations_per_template = 4,
           seed = 42) {
  set.seed(seed) # ensure reproducible randomness

  # grab the transcript location on the genome
  cds <- get_cds(annotation, ensemble_transcript_id)
  mut_genomic <- get_genomic_mutation(cds, mutation_loci)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  aa_cds_seq <- Biostrings::translate(cds_seq)

  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    stop("Your `mutation_original` is not the same as the one on the transcript sequence! Check transcirpt id.")
  }

  genomic_site <- resize(mut_genomic, width = extension * 2, fix = "center")
  genomic_seq <- Biostrings::getSeq(genome, genomic_site)
  names(genomic_seq) <- mutation_name

  # figure out mask - part of the template that actually is not part of the cds!!!
  mask_genome <- intersect(genomic_site, cds[[1]]) # this part on the genome is used for cds
  mask_cds <- ranges(GenomicFeatures::mapToTranscripts(mask_genome, cds)) # this part of the cds we are actually mutating
  # now - make it relative to the genome sequence we operate on
  tx_loci_with_introns <- resize(IRanges(mutation_loci + 1, width = 1),
                                 width = extension * 2, fix = "center")
  mask_seq <- shift(mask_cds, - (start(tx_loci_with_introns) - 1)) # parts of the sequence that are used for the cds

  mutations <- get_all_possilbe_mutations(
    mutation_name, positions_to_mutate, mutation_loci, cds_seq, aa_cds_seq, extension)

  # original mutation
  origin_mutation <- GRanges(
    seqnames = mutation_name,
    ranges = IRanges(extension + 1, width = 1,
                     names = mutation_name),
    strand = "+",
    original = as.character(genomic_seq[[1]][extension + 1]),
    replacement = mutation_replacement, shift = 0,
    codon = ceiling(mutation_loci / 3))

  # find and score guides in a window
  guides <- get_guides_and_scores(origin_mutation, mutation_name, guide_distance, genomic_seq)
  pams <- resize(flank(guides, width = 3, start = FALSE), width = 2, fix = "end") # GG part

  # now the repair template sequences
  repair_template <- GRanges()
  # prepare mutations ranges
  template_range <- promoters(ranges(origin_mutation), upstream = 49, downstream = 51)

  # mutations can't also be using the same codon more than 1 time
  # mutations can't repeat
  # we prefer mutations that mutate the pam > guide > the rest
  for (i in seq_along(guides)) { # each guide gets his own optimized mutations
    muts <- get_combinations_of_mutations_for_guide(mutations, mutations_per_template, pams[i], guides[i])
    mutated_seq <- replaceAt(genomic_seq[[1]],
                             at = ranges(mutations[muts[[1]]]),
                             value = DNAStringSet(mutations[muts[[1]]]$replacement))
    rt <- GRanges(seqnames = mutation_name,
                  ranges = template_range,
                  strand = "+",
                  replacement = as.character(extractAt(mutated_seq, template_range)),
                  pam_disrupted = muts[[2]],
                  guide_disrupted = muts[[3]])
    temp_name <- paste0("Template_", as.character(names(guides[i])),
                        "_Mut_", paste0(muts, collapse = "_"))
    names(rt) <- temp_name
    repair_template <- c(repair_template, rt)
  }


  # first we write the sequence
  Biostrings::writeXStringSet(genomic_seq,
                              file.path(output_dir, paste0(mutation_name, ".fa")))
  # write bed/excel
  origin_mutation$codon <- NULL
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  all_combined <- c(origin_mutation, mutations, guides, repair_template)
  gr <- as.data.frame(all_combined)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  gxt <- data.frame("Guide name" = names(guides),
                    "Guide sequence" = guides$original,
                    "Repair template name " = names(repair_template),
                    "Repair template sequence" = repair_template$replacement,
                    pam_disrupted = repair_template$pam_disrupted,
                    guide_disrupted = repair_template$guide_disrupted,
                    score_rank = guides$rank_by_scores)
  gxt <- gxt[order(gxt$pam_disrupted, gxt$guide_disrupted, -gxt$score_rank, decreasing = T), ]
  write.table(gxt, file = file.path(output_dir, paste0(mutation_name, "_0based_guides_x_templates.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(all_combined, file.path(output_dir, paste0(mutation_name, ".gff3")))
  rtracklayer::export.gff3(GRanges(seqnames = mutation_name,
                                   ranges = mask_seq,
                                   strand = "+",
                                   type = "CDS"), file.path(output_dir, paste0(mutation_name, "_cds.gff3")))
}
