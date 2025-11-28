test_that("process_strand finds guides correctly", {
  mutated_seq <- Biostrings::DNAString(
    "AGCTAGCTAGCTAGAGCTAGCTAGCTAGCTAGCTAGCTAGCTGAGCTAGGAGCTAGCTAGCTAGCTAGCTAGCTAGCTGAGCTATAGCTAGCTAGCTAGCTGAGCTA"
  )
  window_genomic <- GenomicRanges::GRanges("chr1:40-60:+")
  pam_pattern <- "NGG"
  strand_char <- "+"
  scorers <- list()
  score_efficiency <- FALSE

  # Expected result: one guide found
  guides <- process_strand(mutated_seq, pam_pattern, strand_char, scorers, score_efficiency,
                           padding_score = 0, offset_body = 0)
  expect_equal(length(guides), 1)
  expect_equal(as.character(guides$original), "GCTAGCTAGCTAGCTGAGCT")
  expect_equal(as.character(guides$with_pam), "GCTAGCTAGCTAGCTGAGCTAGG")
})
