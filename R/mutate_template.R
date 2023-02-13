# rm(list = ls(all.names = TRUE))
# gc(reset = T)
#
# library(Biostrings)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(IRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)

synonymously_mutate_template <- function(
    ensemble_transcript_id = "ENST00000399837.8",
    mutation_loci = 506,
    mutation_original = "G",
    mutation_replacement = "A",
    mutation_name = "R169Q",
    output_dir = "~/",
    extension = 50,
    positions_to_mutate = -20:20,
    mutations_per_template = 3,
    annotation = "gencode.v42.annotation.gff3",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    n = 96,
    seed = 42)

{
  set.seed(seed) # ensure reproducible randomness

  # grab the transcript location on the genome
  txdb <- GenomicFeatures::makeTxDbFromGFF(annotation)
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = T)
  cds <- cds[!duplicated(names(cds))]
  cds <- cds[ensemble_transcript_id]
  mut_genomic <- GenomicFeatures::pmapFromTranscripts(
    IRanges(mutation_loci, width = 1), cds[ensemble_transcript_id])
  mut_genomic <- mut_genomic[[1]][mut_genomic[[1]]$hit]
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  aa_cds_seq <- Biostrings::translate(cds_seq)

  # # not needed?
  # template <- DNAString("ccccgtccatcagaaaaatgttccaagtggattctgctggaggattatcggaagcgggtgcagaacgtcactgagtttgatgacaggtgagtagtagttc")
  # pairwiseAlignment(genomic_seq, template, type = "overlap")

  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    warning("Your `mutation_original` is not the same as the one on the transcript sequence!")
  }

  # grab sequence +- extension bp around mutation site - might include introns!!!
  genomic_site <- mut_genomic + (extension - 1)
  if (as.character(strand(genomic_site)) == "-") { # on minus we add 1 from the right
    genomic_site <- resize(genomic_site, width = extension * 2, fix = "start")
  } else { # on plus we remove 1 from the right
    genomic_site <- resize(genomic_site, width = extension * 2, fix = "end")
  }
  genomic_seq <- Biostrings::getSeq(genome, genomic_site)
  names(genomic_seq) <- mutation_name

  # 3 extra mutations that are synonymous
  pp <- positions_to_mutate
  # remove positions of our main mutation codon
  this_pos_mutation_loci <- mutation_loci
  codon_count <- ceiling(this_pos_mutation_loci / 3)
  codon_end <- codon_count * 3
  codon_start <- codon_end - 2
  codon_position <- this_pos_mutation_loci - codon_start + 1
  original_codon <- aa_cds_seq[codon_count]
  original_codon_seq <- cds_seq[codon_start:codon_end]

  # prepare mutations ranges
  # original mutation
  if (as.character(strand(genomic_site)) == "-") {
    origin_mutation <- GRanges(
      seqnames = mutation_name,
      ranges = IRanges(extension, width = 1,
                       names = mutation_name),
      strand = "+",
      original = as.character(genomic_seq[[1]][extension + 1]),
      replacement = mutation_replacement, shift = 0,
      codon = codon_count)
  } else {
    origin_mutation <- GRanges(
      seqnames = mutation_name,
      ranges = IRanges(extension + 1, width = 1,
                       names = mutation_name),
      strand = "+",
      original = as.character(genomic_seq[[1]][extension + 1]),
      replacement = mutation_replacement, shift = 0,
      codon = codon_count)
  }


  # 0 is position of the codon_position, therefore
  # filter out other positionsof that codon
  if (codon_position == 3) {
    pp <- pp[!pp %in% -2:0]
  } else if (codon_position == 2) {
    pp <- pp[!pp %in% -1:1]
  } else {
    pp <- pp[!pp %in% 0:2]
  }

  # figure out all possible mutations that are synonymous
  sp <- c()
  mutations <- GRanges()
  for (i in pp) {
    this_pos_mutation_loci <- mutation_loci + i
    codon_count <- ceiling(this_pos_mutation_loci / 3)
    codon_end <- codon_count * 3
    codon_start <- codon_end - 2
    codon_position <- this_pos_mutation_loci - codon_start + 1
    original_codon <- aa_cds_seq[codon_count]
    original_codon_seq <- cds_seq[codon_start:codon_end]

    # we want to change only codon_position and keep same codon
    positions_that_we_keep <- c(1:3)[!c(1:3) %in% codon_position]
    original_codon_seq_ <- original_codon_seq[positions_that_we_keep]
    alternate <- GENETIC_CODE[GENETIC_CODE == as.character(original_codon)]

    # filter out original codon from alternate
    alternate <- alternate[names(alternate) != as.character(original_codon_seq)]

    # filter out those codons that don't have the same extension as original codon
    alt_names <- sapply(names(alternate), function(x) {
      paste0(strsplit(x, "")[[1]][positions_that_we_keep], collapse = "")
    })
    alternate <- alternate[alt_names == as.character(original_codon_seq_)]
    if (length(alternate) == 0) next # original position is crucial to this codon
    # instead of picking randomly we will use all available alternate codons here
    replacement <- sapply(names(alternate), function(x) strsplit(x, "")[[1]][codon_position])
    mut <- GRanges(seqnames = mutation_name,
                   ranges = IRanges(extension + 1 + i, width = 1),
                   strand = "+",
                   original = as.character(original_codon_seq[codon_position]),
                   replacement = "",
                   shift = i,
                   codon = codon_count)
    mut <- rep(mut, length(replacement))
    mut$replacement <- replacement
    mutations <- c(mutations, mut)
  }
  names(mutations) <- seq_along(mutations)

  # from above we need to select 3 mutations that can be unique to each template
  # select 7 groups of 3 mutations
  # each codon can't be reused in this calculation

  # lsit all possible combinations of 3 mutations by index
  tc <- n
  available_muts <- seq_along(mutations)
  selected_combs <- matrix(nrow = tc, ncol = 3)
  i <- 0
  # we just randomize
  stop_counter <- 100000
  stop_counter_i <- 0
  while (i < tc) {
    i <- i + 1
    stop_counter_i <- stop_counter_i + 1
    available_muts <- seq_along(mutations)
    mut123 <- sample(available_muts, size = 3, replace = F)
    mut123 <- sort(mut123)
    # mutations can't also be using the same codon more than 1 time
    # mutations can't repeat
    if ((length(unique(mutations[mut123]$codon)) != 3) |
        (any(apply(selected_combs, 1, function(x) all(x == mut123)), na.rm = T))) {
      i <- i - 1 # this randomization failed, try again
    } else {
      selected_combs[i, ] <- mut123
    }

    if (stop_counter_i == stop_counter) {
      stop("Can't produce that many templates with these settings.")
    }
  }

  # now the repair template sequences
  repair_template <- DNAStringSet()
  for (i in seq_len(tc)) {
    muts <- mutations[selected_combs[i, ]]
    mutated_seq <- replaceAt(genomic_seq[[1]],
                             at = ranges(muts),
                             value = DNAStringSet(muts$replacement))
    temp_name <- paste0("Template_", i, "_Mut_", paste0(selected_combs[i, ], collapse = "_"))
    repair_template[[temp_name]] <- mutated_seq
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
}
