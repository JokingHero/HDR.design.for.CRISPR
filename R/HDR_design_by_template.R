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
                               positions_to_mutate = -19:19,
                               mutations_per_template = 3,
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
  selected_combs <- get_combinations_of_mutations(mutations, 1, mutations_per_template)

  # now the repair template sequences
  repair_template <- GRanges()
  #hdr_probes <- GRanges()
  shifts <- c(0, -30, -20, -10, 10, 20, 30)
  # prepare mutations ranges
  # original mutation
  origin_mutation <- GRanges(
    seqnames = mutation_name,
    ranges = IRanges(extension + 1, width = 1,
                     names = mutation_name),
    strand = "+",
    original = as.character(genomic_seq[[1]][extension + 1]),
    replacement = mutation_replacement, shift = 0,
    codon = ceiling(mutation_loci / 3))
  template_range <- promoters(ranges(origin_mutation), upstream = 49, downstream = 51)
  for (i in 1:7) { # previously we did 7 unique templates as 7 unique mutations, now all 7 get same mutations
    muts <- mutations[selected_combs[1, ]]
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
                        "_Mut_", paste0(selected_combs[1, ], collapse = "_"))
    names(rt) <- temp_name
    repair_template <- c(repair_template, rt)

    # probes <- design_probes(mutation_name, start(origin_mutation) - 18, start(origin_mutation) + 18, mutated_seq,
    #                         len_min = 19, len_max = 25, tmin = 50, tmax = 65,
    #                         origin_mut_start = extension + 1)
    # names(probes) <- paste0("HDR probe ", seq_along(probes), " for ", temp_name)
    # hdr_probes <- c(hdr_probes, probes)
  }

  # find guides in +-100 window CCN/NGG get these sequences
  window <- resize(origin_mutation, width = guide_distance * 2, fix = "center")
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

  # # probes
  # # Reference probes - outside of the SNP area
  # # Length/Tm= 20-25 bp, 60+/-1 oC
  # ref_seq <- IRanges(start = start(origin_mutation), width = 1) + 120
  # ref_seq <- setdiff(ref_seq, IRanges(start = start(origin_mutation), width = 1) + 20)
  # # two chunks that can be used for Ref probes
  # ref_probes <- c(design_probes(mutation_name, start(ref_seq[1]), end(ref_seq[1]), genomic_seq[[1]],
  #                               origin_mut_start = extension + 1),
  #                 design_probes(mutation_name, start(ref_seq[2]), end(ref_seq[2]), genomic_seq[[1]],
  #                               origin_mut_start = extension + 1))
  # names(ref_probes) <- paste0("Ref probe ", seq_along(ref_probes))
  #
  # # NHEJ probe - ~20 bp and have Tm values 58-60oC
  # # contain the original mutated sequence  no SNPs and no correction of the mutation
  # nhej_probes <- design_probes(mutation_name,
  #                              start(origin_mutation) - 18, start(origin_mutation) + 18, mutated_seq,
  #                              len_min = 19, len_max = 25, tmin = 50, tmax = 65,
  #                              origin_mut_start = extension + 1)
  # names(nhej_probes) <- paste0("NHEJ probe ", seq_along(nhej_probes))


  # first we write the sequence
  Biostrings::writeXStringSet(genomic_seq,
                              file.path(output_dir, paste0(mutation_name, ".fa")))
  # write bed/excel
  origin_mutation$codon <- NULL
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  all_combined <- c(origin_mutation, mutations, guides, repair_template)
                    #hdr_probes, nhej_probes, ref_probes)
  gr <- as.data.frame(all_combined)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  gxt <- data.frame("Guide name" = rep(names(guides), each = length(repair_template)),
                    "Guide sequence" = rep(guides$original, each = length(repair_template)),
                    "Repair template name " = rep(names(repair_template), length(guides)),
                    "Repair template sequence" = rep(repair_template$replacement, length(guides)))
  write.table(gxt, file = file.path(output_dir, paste0(mutation_name, "_0based_guides_x_templates.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(all_combined, file.path(output_dir, paste0(mutation_name, ".gff3")))
  rtracklayer::export.gff3(GRanges(seqnames = mutation_name,
                                   ranges = mask_seq,
                                   strand = "+",
                                   type = "CDS"), file.path(output_dir, paste0(mutation_name, "_cds.gff3")))
}
