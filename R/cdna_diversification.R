#' @title Get genomic seqeunce as a fasta file
#'
#' @description Outputs genomic sequence and homology arms of the specified lengths
#'
#' @param chromosome What is the chromosome of your target site?
#' @param position Position on the chromosome
#' @param bp How many bp of the seqeunce around the cut site you would like to extract?
#' @param homology_arm_length How long should homology arms be? Right arm includes `position`.
#' @param output_file To which file output the results? Writes fasta file. By default it does not save the file.
#' @param genome Which genome to use? By default it uses BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38.
#' @return returns the DNAStringSet of the extracted sequences
#' @import Biostrings GenomicRanges IRanges BSgenome.Hsapiens.UCSC.hg38
#' @export
#'
get_genomic_seq <- function(
    chromosome = "chr1",
    position = 2313,
    bp = 3000,
    homology_arm_length = 400,
    output_file = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) {

  cut_site <- GenomicRanges::GRanges(seqnames = chromosome,
                                     ranges = IRanges::IRanges(
                                       start = position, width = 1),
                                     strand = "+")
  GenomeInfoDb::seqlevelsStyle(cut_site) <- GenomeInfoDb::seqlevelsStyle(genome)
  GenomeInfoDb::seqlevels(cut_site) <- GenomeInfoDb::seqlevels(genome)
  GenomeInfoDb::seqlengths(cut_site) <- GenomeInfoDb::seqlengths(genome)
  cut_site_extended <- resize(cut_site, width = bp * 2, fix = "center")
  cut_site_extended <- trim(cut_site_extended)
  cut_site_extended_seq <- Biostrings::getSeq(genome, cut_site_extended)

  homology_arm_left <- flank(cut_site, width = homology_arm_length, start = T)
  homology_arm_left <- trim(homology_arm_left)
  homology_arm_left_seq <- Biostrings::getSeq(genome, homology_arm_left)

  homology_arm_right <- flank(cut_site, width = homology_arm_length, start = F)
  homology_arm_right <- trim(homology_arm_right)
  homology_arm_right_seq <- Biostrings::getSeq(genome, homology_arm_right)

  seq <- c(cut_site_extended_seq, homology_arm_left_seq, homology_arm_right_seq)
  names(seq) <- c(paste0(cut_site_extended),
                  paste0(homology_arm_left),
                  paste0(homology_arm_right))
  if (!is.null(output_file)) {
    writeXStringSet(seq, filepath = output_file, compress = F)
  }
  return(seq)
}

#' @title Design templates with innocent synonymous SNPs
#'
#' @description Design templates with mutations that will not change coding sequence
#'
#' @param ensemble_transcript_id This has to be ensemble transcript id e.g. ENST00000307851.9
#' @param annotation File path to the annotation file, a gff3.
#' @param output_file Path and name of the out fasta file e.g "~/mutated_cds.fa"
#' @param genome BSgenome of your genome, compatible with annotation file, by default it is hg38.
#' @param seed Ensures reproducibility of the random parts of the pipeline.
#' @return writes file to `output_file`, might overwrite
#' @import Biostrings GenomicFeatures GenomicRanges IRanges BSgenome.Hsapiens.UCSC.hg38
#' @export
#'
synonymously_mutate_cds <- function(
    ensemble_transcript_id,
    annotation,
    output_file,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    seed = 42) {

  set.seed(seed) # ensure reproducible randomness

  # grab the transcript location on the genome
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(annotation))
  cds <- get_cds(txdb, ensemble_transcript_id)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  aa_cds_seq <- Biostrings::translate(cds_seq)

  codon_count <- length(aa_cds_seq)

  # figure out all possible mutations that are synonymous
  synonymous_cds <- c()
  for (i in seq_len(codon_count)) {
    alternate <- Biostrings::GENETIC_CODE[Biostrings::GENETIC_CODE == as.character(aa_cds_seq[i])]

    # filter out original codon from alternate
    original_cds_seq <- as.character(cds_seq[(i*3 - 2):(i*3)])
    alternate <- alternate[names(alternate) != original_cds_seq]

    # is this needed?!
    # # filter out those codons that don't have the same extension as original codon
    # alt_names <- sapply(names(alternate), function(x) {
    #   paste0(strsplit(x, "")[[1]][positions_that_we_keep], collapse = "")
    # })
    # alternate <- alternate[alt_names == as.character(original_codon_seq_)]
    if (length(alternate) == 0) {
      synonymous_cds <- c(synonymous_cds, original_cds_seq)
    } else {
      synonymous_cds <- c(synonymous_cds, sample(names(alternate), 1))
    }
  }

  synonymous_cds <- DNAString(paste0(synonymous_cds, collapse = ""))

  aa_synonymous_cds <- Biostrings::translate(synonymous_cds)
  # verify the same that synonymous changes were introduced
  if (aa_synonymous_cds != aa_cds_seq) {
    stop("Not all sequences seem to be perfectly non-synonymous changes. This is a bug, report to the authors.")
  }

  if (as.character(aa_synonymous_cds[codon_count]) != "*") {
    aa_synonymous_cds <- paste0(c(as.character(synonymous_cds), "TAG"), collapse = "")
  }

  synonymous_cds <- DNAStringSet(synonymous_cds)
  names(synonymous_cds) <- paste0(ensemble_transcript_id, " synonymous mutated with seed ", seed)

  # first we write the sequence
  Biostrings::writeXStringSet(synonymous_cds, output_file)
}
