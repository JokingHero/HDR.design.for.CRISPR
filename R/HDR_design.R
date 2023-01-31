# rm(list = ls(all.names = TRUE))
# gc(reset = T)
#
# set.seed(42) # ensure reproducible randomness
# library(Biostrings)
# library(GenomicFeatures)
# library("BSgenome.Hsapiens.UCSC.hg38")
# genome <- BSgenome.Hsapiens.UCSC.hg38
#
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
#
# fasta_output_file <- paste0(mutation_name, ".fa")
# coordinates_output_file <- paste0(mutation_name, ".bed")
# extension <- 400
#
# # grab the transcript location on the genome
# txdb <- GenomicFeatures::makeTxDbFromGFF("gencode.v42.annotation.gff3")
# cds <- cdsBy(txdb, by = "tx", use.names = T)
# cds <- cds[!duplicated(names(cds))]
# cds <- cds[ensemble_transcript_id]
# mut_genomic <- pmapFromTranscripts(IRanges(mutation_loci, width = 1), cds[ensemble_transcript_id])
# mut_genomic <- mut_genomic[[1]][mut_genomic[[1]]$hit]
# cds_seq <- extractTranscriptSeqs(genome, cds)[[1]]
# aa_cds_seq <- translate(cds_seq)
#
# if (as.character(getSeq(genome, mut_genomic)) != mutation_original) {
#   warning("Your `mutation_original` is not the same as the one on the transcript sequence!")
# }
#
# # grab sequence +- 400bp around mutation site
# genomic_site <- mut_genomic + extension
# genomic_seq <- getSeq(genome, genomic_site)
# names(genomic_seq) <- mutation_name
#
# # prepare mutations ranges
# # original mutation
# mutations <- GRanges(seqnames = mutation_name,
#         ranges = IRanges(extension + 1, width = 1,
#                          names = mutation_name),
#         strand = "+",
#         original = as.character(genomic_seq[[1]][extension + 1]),
#         replacement = mutation_replacement, shift = 0)
#
# # 3 extra mutations that are synonymous
# pp <- -10:10
# # remove positions of our main mutation codon
# this_pos_mutation_loci <- mutation_loci
# codon_count <- ceiling(this_pos_mutation_loci / 3)
# codon_end <- codon_count * 3
# codon_start <- codon_end - 2
# codon_position <- this_pos_mutation_loci - codon_start + 1
# original_codon <- aa_cds_seq[codon_count]
# original_codon_seq <- cds_seq[codon_start:codon_end]
# # 0 is position of the codon_position, therefore
# if (codon_position == 3) {
#   pp <- pp[!pp %in% -2:0]
# } else if (codon_position == 2) {
#   pp <- pp[!pp %in% -1:1]
# } else {
#   pp <- pp[!pp %in% 0:2]
# }
#
# # randomize order for fun
# pp <- sample(pp, length(pp), replace = F)
# sp <- c()
# locked_pos <- c()
# for (i in pp) {
#   if (i %in% locked_pos) next
#   this_pos_mutation_loci <- mutation_loci + i
#   codon_count <- ceiling(this_pos_mutation_loci / 3)
#   codon_end <- codon_count * 3
#   codon_start <- codon_end - 2
#   codon_position <- this_pos_mutation_loci - codon_start + 1
#   original_codon <- aa_cds_seq[codon_count]
#   original_codon_seq <- cds_seq[codon_start:codon_end]
#
#   # we want to change only codon_position and keep same codon
#   positions_that_we_keep <- c(1:3)[!c(1:3) %in% codon_position]
#   original_codon_seq_ <- original_codon_seq[positions_that_we_keep]
#   alternate <- GENETIC_CODE[GENETIC_CODE == as.character(original_codon)]
#
#   # filter out original codon from alternate
#   alternate <- alternate[names(alternate) != as.character(original_codon_seq)]
#
#   # filter out those codons that don't have the same extension as original codon
#   alt_names <- sapply(names(alternate), function(x) {
#     paste0(strsplit(x, "")[[1]][positions_that_we_keep], collapse = "")
#   })
#   alternate <- alternate[alt_names == as.character(original_codon_seq_)]
#   if (length(alternate) == 0) next # original position is crucial to this codon
#   # now we pick randomly
#   syn_codon <- sample(names(alternate), size = 1)
#
#   mutations <- c(
#     mutations,
#     GRanges(seqnames = mutation_name,
#             ranges = IRanges(extension + 1 + i, width = 1,
#                              names = paste0("Mutation_", as.character(i))),
#             strand = "+",
#             original = as.character(original_codon_seq[codon_position]),
#             replacement = strsplit(syn_codon, "")[[1]][codon_position],
#             shift = i))
#   # i is codon_position
#   if (codon_position == 1) {
#     lock <- c(i, i + 1, i + 2)
#   } else if (codon_position == 2) {
#     lock <- c(i - 1, i, i + 2)
#   } else {
#     lock <- c(i - 2, i - 1, i)
#   }
#   locked_pos <- c(locked_pos, lock)
# }
#
#
# # find guides in +-100 window CCN/NGG get these sequences
# window <- mutations[1] + 50
# pam_fwd <- IRanges(matchPattern(DNAString("GG"), genomic_seq[[1]]))
# pam_rve <- IRanges(matchPattern(DNAString("CC"), genomic_seq[[1]]))
# # restrict to window
# pam_fwd <- pam_fwd[pam_fwd %over% ranges(window)]
# pam_rve <- pam_rve[pam_rve %over% ranges(window)]
#
# # now we make just guides
# pam_fwd <- flank(pam_fwd + 1, width = 20, start = T)
# names(pam_fwd) <- paste0("NGG_", as.character(seq_along(pam_fwd)))
# pam_fwd <- GRanges(seqnames = mutation_name,
#                    ranges = pam_fwd,
#                    strand = "+",
#                    original = as.character(extractAt(genomic_seq[[1]], pam_fwd)),
#                    replacement = "",
#                    shift = end(pam_fwd) - start(mutations[1]))
#
# pam_rve <- flank(pam_rve + 1, width = 20, start = F)
# names(pam_rve) <- paste0("CCN_", as.character(seq_along(pam_rve)))
# pam_rve <- GRanges(seqnames = mutation_name,
#                    ranges = pam_rve,
#                    strand = "-",
#                    original = as.character(reverseComplement(extractAt(genomic_seq[[1]], pam_rve))),
#                    replacement = "",
#                    shift = start(pam_rve) - start(mutations[1]))
# guides <- c(pam_fwd, pam_rve)
#
#
# # now the repair template sequences
# # we restrict to random 3 mutations
# mutations <- mutations[c(1, sample(2:length(mutations), 3))]
# not_first <- 2:length(mutations)
# mutated_seq <- replaceAt(genomic_seq[[1]],
#                          at = ranges(mutations[not_first]),
#                          value = DNAStringSet(mutations[not_first]$replacement))
# # now the repair templates
# repair_template <- promoters(ranges(mutations[1]), upstream = 49, downstream = 51)
# shifts <- c(0, -30, -20, -10, 10, 20, 30)
# repair_template <- c(repair_template,
#                      shift(repair_template, shift = -30),
#                      shift(repair_template, shift = -20),
#                      shift(repair_template, shift = -10),
#                      shift(repair_template, shift = 10),
#                      shift(repair_template, shift = 20),
#                      shift(repair_template, shift = 30))
# names(repair_template) <- paste0("Template_", as.character(shifts))
# repair_template <-
#   GRanges(seqnames = mutation_name,
#         ranges = repair_template,
#         strand = "+",
#         original = as.character(extractAt(genomic_seq[[1]], repair_template)),
#         replacement = as.character(extractAt(mutated_seq, repair_template)),
#         shift = shifts)
#
# writeXStringSet(genomic_seq, fasta_output_file)
# # write bed/excel
# # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
# gr <- as.data.frame(c(mutations, guides, repair_template))
# gr$start <- gr$start - 1
# write.table(gr, file = coordinates_output_file,
#             quote = F, sep = "\t", row.names = F, col.names = T)
