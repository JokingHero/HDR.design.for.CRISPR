# Tests for functions in R/design_hdr.R
annot <- system.file(
  "data", "gencode.v42.annotation_mini.gff3", package = "HDR.design.for.CRISPR")
txdb <- suppressMessages(
  suppressWarnings(txdbmaker::makeTxDbFromGFF(annot, format = "gff3")))

test_that("create_template_and_probes works as expected", {
  # Mock inputs
  muts <- list(
    mutations = GRanges("chr1:15-15", REF = "C", ALT = "T"),
    pam_disrupted_count = 1,
    guide_disrupted_count = 0,
    total_cadd = 10,
    any_overlaps_noncoding = FALSE,
    total_compatibility_score = 5
  )
  guide_name <- "TestGuide"
  template_range <- GRanges("chr1:1-30")
  template_ref <- DNAString(paste(rep("A", 30), collapse = ""))
  design_name <- "TestDesign"
  do_probes <- FALSE
  probe_params <- list()

  # Call the function
  result <- create_template_and_probes(
    muts, guide_name, template_range, template_ref, design_name, do_probes, probe_params
  )

  # Assertions
  expect_true("template" %in% names(result))
  expect_true("probes" %in% names(result))
  expect_equal(length(result$template), 1)
  expect_equal(result$template$pam_disrupted, 1)
  expect_true(grepl("Template_TestGuide_Mut_", names(result$template)))
  # Since do_probes is FALSE, expect empty probes
  expect_equal(length(result$probes), 0)
})
