# rm(list = ls(all.names = TRUE))
# gc(reset = T)
#
# library(Biostrings)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(IRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
#
# ensemble_transcript_id = "ENST00000399837.8"
# mutation_loci = 506
# mutation_original = "G"
# mutation_replacement = "A"
# mutation_name = "R169Q"
# output_dir = "~/"
# annotation = "/mnt/corsair/Projects/uib/CRIPSR/Oslo/Oline_seqeuncing_HDR_14_07_2023/gencode.v42.annotation.gff3.gz"
# genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# extension = 50
# positions_to_mutate = -30:30
# mutations_per_template = 4
# n = 96
# seed = 42

#' @title Design templates with innocent synonymous SNPs
#'
#' @description Design templates with mutations that will not change coding sequence
#'
#' @param ensemble_transcript_id This has to be ensemble transcript id e.g. ENST00000307851.9
#' @param mutation_loci Position on the gene that is mutated e.g. 245
#' @param mutation_original What was the original base on the genome? e.g. A
#' @param mutation_replacement What is the mutated base? e.g. G
#' @param mutation_name What is your name for this mutation e.g. tyr82cys
#' @param output_dir Where the files should be generated? Make sure you have writing permissions there.
#' @param annotation File path to the annotation file, a gff3.
#' @param genome BSgenome of your genome, compatible with annotation file, by default it is hg38.
#' @param extension How many bases upstream/downstream from the mutation loci should we include Default is 50.
#' @param positions_to_mutate Which positions from the mutation are available for mutation. By default its -30:30
#' leaving 20bp on each side of the template for the in case of incomplete integration of the template. Also codon
#' occupied by the mutation is not included in the change.
#' @param mutations_per_template How many codons should be mutated per template.
#' @param n How many templates to generate? Default is set to 96.
#' @param seed Ensures reproducibility of the random parts of the pipeline.
#' @return writes files to the specified directory, might overwrite
#' @import Biostrings GenomicFeatures GenomicRanges IRanges BSgenome.Hsapiens.UCSC.hg38
#' @importFrom utils write.table
#' @export
#'
synonymously_mutate_template <- function(
    ensemble_transcript_id,
    mutation_loci,
    mutation_original,
    mutation_replacement,
    mutation_name,
    output_dir,
    annotation,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    extension = 50,
    positions_to_mutate = -30:30,
    mutations_per_template = 3,
    n = 96,
    seed = 42) {

  set.seed(seed) # ensure reproducible randomness

  # grab the transcript location on the genome
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(annotation))
  cds <- get_cds(txdb, ensemble_transcript_id)
  mut_genomic <- get_genomic_mutation(cds, mutation_loci)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  aa_cds_seq <- Biostrings::translate(cds_seq)

  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    stop("Your `mutation_original` is not the same as the one on the transcript sequence! Check transcirpt id.")
  }

  # grab sequence +- extension bp around mutation site - might include introns!!!
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
  selected_combs <- get_combinations_of_mutations(mutations, n, mutations_per_template)

  # now the repair template sequences
  repair_template <- DNAStringSet()
  for (i in seq_len(nrow(selected_combs))) {
    muts <- mutations[selected_combs[i, ]]
    seq_to_mut <- genomic_seq[[1]]
    mutated_seq <- replaceAt(seq_to_mut,
                             at = ranges(muts),
                             value = DNAStringSet(muts$replacement))
    temp_name <- paste0("Template_", i, "_Mut_", paste0(selected_combs[i, ], collapse = "_"))
    repair_template[[temp_name]] <- mutated_seq
  }

  # lets do some verification that all codons work
  if (!all(sapply(repair_template, function(x) {
    translate(replaceAt(cds_seq, at = mask_cds, value = x[mask_seq])) == aa_cds_seq
  }))) {
    stop("Not all sequences seem to be perfectly non-synonymous changes. This is a bug, report to authors.")
  }

  # first we write the sequence
  Biostrings::writeXStringSet(repair_template,
                              file.path(output_dir, paste0(mutation_name, "_templates.fa")))
  # write bed/excel
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  gr <- as.data.frame(mutations)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based_mutations.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(mutations, file.path(output_dir, paste0(mutation_name, ".gff3")))
  rtracklayer::export.gff3(GRanges(seqnames = mutation_name,
                                   ranges = mask_seq,
                                   strand = "+",
                                   type = "CDS"), file.path(output_dir, paste0(mutation_name, "_cds.gff3")))
}
