# Tests for functions in R/design_hdr.R
annot <- testthat::test_path("testdata", "gencode.v42.annotation_mini.gff3")

# library(HDR.design.for.CRISPR)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# library(GenomicScores)

test_that("End-to-end test for STAT3 (minus strand)", {
  expected_files_dir <- testthat::test_path("testdata", "STAT3")
  output_dir <- tempfile(pattern = "STAT3_minus_strand_test_")
  dir.create(output_dir)

  # Using a known transcript and mutation for STAT3
  design_hdr(
    design_name = "STAT3",
    optimization_scheme = "balanced",
    seed = 42,
    chrom = "chr17",
    variant_start = 42346635,
    variant_end = 42346635,
    REF = "G",
    ALT = "A",
    ALT_on_guides = TRUE,
    ALT_on_templates = FALSE,
    output_dir = output_dir,
    annotation = annot,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    maximum_mutations_per_template = 3,
    snps = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
    cadd = GenomicScores::getGScores("cadd.v1.6.hg38")
  )

  # Get lists of files in both directories
  # we only grab CSV because they don't save dates
  output_files <- list.files(output_dir, full.names = TRUE, pattern = "*.csv")
  expected_files <- list.files(expected_files_dir, full.names = TRUE, pattern = "*.csv")

  # Extract just the filenames for comparison
  output_basenames <- basename(output_files)
  expected_basenames <- basename(expected_files)

  # 1. Check if the set of filenames is the same
  expect_equal(sort(output_basenames), sort(expected_basenames),
    info = "Filenames in output directory should match expected filenames."
  )

  # 2. Compare content of each file
  if (length(output_basenames) > 0) {
    for (filename in output_basenames) {
      filename <- output_basenames[7]
      output_filepath <- file.path(output_dir, filename)
      expected_filepath <- file.path(expected_files_dir, filename)

      # Ensure the expected file exists
      expect_true(file.exists(expected_filepath),
        info = paste("Expected file", filename, "does not exist.")
      )

      output_content <- readLines(output_filepath, warn = FALSE)
      expected_content <- readLines(expected_filepath, warn = FALSE)

      expect_equal(output_content, expected_content,
        info = paste("Content of file", filename, "should match expected content.")
      )
    }
  } else {
    warning("No files were generated in the output directory to compare.")
    expect_gt(length(output_files), 0)
  }

  # Clean up the temporary directory
  unlink(output_dir, recursive = TRUE)
})

test_that("End-to-end test for CTNNB1 (plus strand)", {
  expected_files_dir <- testthat::test_path("testdata", "CTNNB1")
  output_dir <- tempfile(pattern = "CTNNB1_plus_strand_test_")
  dir.create(output_dir)

  # Using a known transcript and mutation for CTNNB1
  # The gencode.v42.annotation_mini.gff3 file confirms that CTNNB1 is on the plus strand.
  design_hdr(
    design_name = "CTNNB1",
    chrom = "chr3",
    variant_start = 41224622,
    variant_end = 41224622,
    REF = "C",
    ALT = "T",
    ALT_on_guides = TRUE,
    ALT_on_templates = FALSE,
    output_dir = output_dir,
    annotation = annot,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    maximum_mutations_per_template = 3,
    snps = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
    cadd = GenomicScores::getGScores("cadd.v1.6.hg38")
  )

  # Get lists of files in both directories
  output_files <- list.files(output_dir, full.names = TRUE, pattern = "*.csv")
  expected_files <- list.files(expected_files_dir, full.names = TRUE, pattern = "*.csv")

  # Extract just the filenames for comparison
  output_basenames <- basename(output_files)
  expected_basenames <- basename(expected_files)

  # 1. Check if the set of filenames is the same
  expect_equal(sort(output_basenames), sort(expected_basenames),
    info = "Filenames in output directory should match expected filenames."
  )

  # 2. Compare content of each file
  if (length(output_basenames) > 0) {
    for (filename in output_basenames) {
      output_filepath <- file.path(output_dir, filename)
      expected_filepath <- file.path(expected_files_dir, filename)

      # Ensure the expected file exists
      expect_true(file.exists(expected_filepath),
        info = paste("Expected file", filename, "does not exist.")
      )

      # Read and compare file contents
      # Assuming text files. For binary files, you might need a different approach.
      output_content <- readLines(output_filepath, warn = FALSE)
      expected_content <- readLines(expected_filepath, warn = FALSE)

      expect_equal(output_content, expected_content,
        info = paste("Content of file", filename, "should match expected content.")
      )
    }
  } else {
    warning("No files were generated in the output directory to compare.")
    expect_gt(length(output_files), 0)
  }

  # Clean up the temporary directory
  unlink(output_dir, recursive = TRUE)
})
