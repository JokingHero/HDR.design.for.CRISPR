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

#' @title Design HDR probes and CRISPR guides
#'
#' @description Design HDR probes, CRISPR guides and NHEJ probes for a given
#' template transcript.
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
#' @param extension How many bases upstream/downstream from the mutation loci should we include Default is 400.
#' @param seed Ensures reproducibility of the random parts of the pipeline.
#' @return writes files to the specified directory, might overwrite
#' @import Biostrings GenomicFeatures GenomicRanges IRanges BSgenome.Hsapiens.UCSC.hg38
#' @importFrom utils write.table
#' @export
#'
design_by_template <- function(ensemble_transcript_id,
                               mutation_loci,
                               mutation_original,
                               mutation_replacement,
                               mutation_name,
                               output_dir,
                               annotation,
                               genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                               guide_distance = 10 + 17,
                               extension = 400,
                               seed = 42) {
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

  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    warning("Your `mutation_original` is not the same as the one on the transcript sequence!")
  }

  # grab sequence +- 400bp around mutation site
  genomic_site <- mut_genomic + extension
  genomic_seq <- Biostrings::getSeq(genome, genomic_site)
  names(genomic_seq) <- mutation_name

  # 3 extra mutations that are synonymous
  pp <- -19:19
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
  origin_mutation <- GRanges(
    seqnames = mutation_name,
    ranges = IRanges(extension + 1, width = 1,
                     names = mutation_name),
    strand = "+",
    original = as.character(genomic_seq[[1]][extension + 1]),
    replacement = mutation_replacement, shift = 0,
    codon = codon_count)

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
  tc <- 7
  available_muts <- seq_along(mutations)
  selected_combs <- matrix(nrow = tc, ncol = 3)
  i <- 1
  while(length(available_muts) > 2 & i <= tc) {
    # in each iteration make the best selection, for harder and harder problem
    # first, the best starting codon will be the one with highest uniqness
    codons_by_counts <- table(mutations[available_muts]$codon)
    codon1 <- names(which.min(codons_by_counts))
    # now randomize selection for that codon
    mut1 <- available_muts[mutations[available_muts]$codon == codon1][1]
    # remove this selection
    available_muts <- available_muts[available_muts != mut1]

    # as we have secured one truly unique codon for our templates
    # we can select one of the less unique codons
    codons_by_counts <- table(mutations[available_muts]$codon)
    codons_by_counts <- codons_by_counts[names(codons_by_counts) != codon1]
    codon2 <- names(which.max(codons_by_counts))
    # now randomize selection for that codon
    mut2 <- available_muts[mutations[available_muts]$codon == codon2][1]
    # remove this selection
    available_muts <- available_muts[available_muts != mut2]

    # now we want to select codon that has big space but not the same as the codon2
    codons_by_counts <- table(mutations[available_muts]$codon)
    codons_by_counts <- codons_by_counts[!names(codons_by_counts) %in% c(codon1, codon2)]
    codon3 <- names(which.max(codons_by_counts))
    # now randomize selection for that codon
    mut3 <- available_muts[mutations[available_muts]$codon == codon3][1]
    # remove this selection
    available_muts <- available_muts[available_muts != mut3]
    selected_combs[i, ] <- c(mut1, mut2, mut3)
    i = i + 1
  }

  # now we fill up the rest with leftovers
  if ((i < tc + 1) & length(available_muts) > 0) {
    mut1 <- available_muts[1]
    available_muts <- available_muts[available_muts != mut1]
    if (length(available_muts) > 0) {
      mut2 <- available_muts[1]
      mut3 <- selected_combs[1, 1] # take the path of least resistance

    } else {
      mut2 <- selected_combs[1, 1]
      mut3 <- selected_combs[2, 1]
    }
    selected_combs[i, ] <- c(mut1, mut2, mut3)
  }

  # now we just randomize
  selected_combs <- Rfast::sort_mat(selected_combs, by.row = T, descending = F)
  while (i < tc) {
    i <- i + 1
    available_muts <- seq_along(mutations)
    mut123 <- sample(available_muts, size = 3, replace = F)
    mut123 <- sort(mut123)
    if (any(apply(selected_combs, 1,function(x) all(x == mut123)), na.rm = T)) {
      i <- i - 1 # this randomization failed, try again
    } else {
      selected_combs[i, ] <- mut123
    }
  }

  # now the repair template sequences
  repair_template <- GRanges()
  hdr_probes <- GRanges()
  shifts <- c(0, -30, -20, -10, 10, 20, 30)
  template_range <- promoters(ranges(origin_mutation), upstream = 49, downstream = 51)
  for (i in seq_len(tc)) {
    muts <- mutations[selected_combs[i, ]]
    mutated_seq <- replaceAt(genomic_seq[[1]],
                             at = ranges(muts),
                             value = DNAStringSet(muts$replacement))
    tri <- shift(template_range, shift = shifts[i])
    rt <- GRanges(seqnames = mutation_name,
                  ranges = tri,
                  strand = "+",
                  replacement = as.character(extractAt(mutated_seq, tri)),
                  shift = shifts[i])
    temp_name <- paste0("Template_", as.character(shifts[i]),
                        "_Mut_", paste0(selected_combs[i, ], collapse = "_"))
    names(rt) <- temp_name
    repair_template <- c(repair_template, rt)

    probes <- design_probes(mutation_name, start(origin_mutation) - 18, start(origin_mutation) + 18, mutated_seq,
                            len_min = 19, len_max = 25, tmin = 50, tmax = 65,
                            origin_mut_start = extension + 1)
    names(probes) <- paste0("HDR probe ", seq_along(probes), " for ", temp_name)
    hdr_probes <- c(hdr_probes, probes)
  }

  # find guides in +-100 window CCN/NGG get these sequences
  window <- origin_mutation + guide_distance
  mutated_seq <- replaceAt(genomic_seq[[1]],
                           at = ranges(origin_mutation),
                           value = DNAStringSet(origin_mutation$replacement))
  pam_fwd <- IRanges(matchPattern(DNAString("GG"), mutated_seq))
  pam_rve <- IRanges(matchPattern(DNAString("CC"), mutated_seq))
  # restrict to window
  pam_fwd <- pam_fwd[pam_fwd %over% ranges(window)]
  pam_rve <- pam_rve[pam_rve %over% ranges(window)]

  # now we make just guides
  pam_fwd <- flank(pam_fwd + 1, width = 20, start = T)
  if (length(pam_fwd) > 0) {
    names(pam_fwd) <- paste0("NGG_", as.character(seq_along(pam_fwd)))
    pam_fwd <- GRanges(seqnames = mutation_name,
                       ranges = pam_fwd,
                       strand = "+",
                       original = as.character(extractAt(mutated_seq, pam_fwd)),
                       shift = end(pam_fwd) - start(origin_mutation))
  }
  pam_rve <- flank(pam_rve + 1, width = 20, start = F)
  if (length(pam_rve) > 0) {
    names(pam_rve) <- paste0("CCN_", as.character(seq_along(pam_rve)))
    pam_rve <- GRanges(seqnames = mutation_name,
                       ranges = pam_rve,
                       strand = "-",
                       original = as.character(reverseComplement(extractAt(mutated_seq, pam_rve))),
                       shift = start(pam_rve) - start(origin_mutation))
  }
  guides <- if (length(pam_fwd) > 0 & length(pam_rve) > 0) {
    c(pam_fwd, pam_rve)
  } else if (length(pam_fwd) > 0) {
    pam_fwd
  } else {
    pam_rve
  }

  # probes
  # Reference probes - outside of the SNP area
  # Length/Tm= 20-25 bp, 60+/-1 oC
  ref_seq <- IRanges(start = start(origin_mutation), width = 1) + 120
  ref_seq <- setdiff(ref_seq, IRanges(start = start(origin_mutation), width = 1) + 20)
  # two chunks that can be used for Ref probes
  ref_probes <- c(design_probes(mutation_name, start(ref_seq[1]), end(ref_seq[1]), genomic_seq[[1]],
                                origin_mut_start = extension + 1),
                  design_probes(mutation_name, start(ref_seq[2]), end(ref_seq[2]), genomic_seq[[1]],
                                origin_mut_start = extension + 1))
  names(ref_probes) <- paste0("Ref probe ", seq_along(ref_probes))

  # NHEJ probe - ~20 bp and have Tm values 58-60oC
  # contain the original mutated sequence  no SNPs and no correction of the mutation
  nhej_probes <- design_probes(mutation_name,
                               start(origin_mutation) - 18, start(origin_mutation) + 18, mutated_seq,
                               len_min = 19, len_max = 25, tmin = 50, tmax = 65,
                               origin_mut_start = extension + 1)
  names(nhej_probes) <- paste0("NHEJ probe ", seq_along(nhej_probes))


  # first we write the sequence
  Biostrings::writeXStringSet(genomic_seq,
                              file.path(output_dir, paste0(mutation_name, ".fa")))
  # write bed/excel
  origin_mutation$codon <- NULL
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  all_combined <- c(origin_mutation, mutations, guides, repair_template,
                    hdr_probes, nhej_probes, ref_probes)
  gr <- as.data.frame(all_combined)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(all_combined, file.path(output_dir, paste0(mutation_name, ".gff3")))
}
