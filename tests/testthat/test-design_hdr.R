# Tests for functions in R/design_hdr.R

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

test_that("design_hdr runs with 'optimal_per_guide' strategy", {
  # Mocking dependencies
  # We assume these functions are tested elsewhere and return expected objects
  mockery::stub(design_hdr, 'prepare_candidate_snps', data.frame())
  mockery::stub(design_hdr, 'get_guides_and_scores_refactored', GRanges("chr1:1-1"))
  mockery::stub(design_hdr, 'find_mutation_combinations', list(mutations = GRanges()))
  mockery::stub(design_hdr, 'create_template_and_probes', list(template = GRanges(), probes = GRanges()))
  mockery::stub(design_hdr, 'export_design_results', NULL) # Don't write files
  mockery::stub(design_hdr, 'getSeq', Biostrings::DNAString("A"))
  mockery::stub(design_hdr, 'txdbmaker::makeTxDbFromGFF', suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(GRanges())))


  # Create a temporary directory for output
  temp_dir <- tempdir()

  # Minimal parameters to run the function
  # Most of these are placeholders as we are mocking the core logic
  expect_invisible(design_hdr(
    design_name = "test_run",
    chrom = "chr1",
    variant_start = 100,
    variant_end = 100,
    REF = "A",
    ALT = "G",
    ALT_on_guides = FALSE,
    ALT_on_templates = TRUE,
    output_dir = temp_dir,
    annotation = "mock.gff", # Mocked to not matter
    strategy = "optimal_per_guide",
    do_probes = FALSE
  ))
})

test_that("design_hdr runs with 'optimal_for_all' and 'all_per_guide' strategies", {
  # Mocking dependencies
  mockery::stub(design_hdr, 'prepare_candidate_snps', data.frame())
  mockery::stub(design_hdr, 'get_guides_and_scores_refactored', GRanges("chr1:1-1"))
  mockery::stub(design_hdr, 'find_mutation_combinations', list(mutations = GRanges()))
  mockery::stub(design_hdr, 'create_template_and_probes', list(template = GRanges(), probes = GRanges()))
  mockery::stub(design_hdr, 'export_design_results', NULL)
  mockery::stub(design_hdr, 'getSeq', Biostrings::DNAString("A"))
  mockery::stub(design_hdr, 'txdbmaker::makeTxDbFromGFF', suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(GRanges())))

  temp_dir <- tempdir()
  common_args <- list(
    design_name = "test_run",
    chrom = "chr1",
    variant_start = 100,
    variant_end = 100,
    REF = "A",
    ALT = "G",
    ALT_on_guides = FALSE,
    ALT_on_templates = TRUE,
    output_dir = temp_dir,
    annotation = "mock.gff",
    do_probes = FALSE
  )

  # Test 'optimal_for_all'
  args_optimal_all <- c(common_args, strategy = "optimal_for_all")
  expect_invisible(do.call(design_hdr, args_optimal_all))

  # Test 'all_per_guide'
  args_all_per_guide <- c(common_args, strategy = "all_per_guide")
  expect_invisible(do.call(design_hdr, args_all_per_guide))
})