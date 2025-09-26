test_that("write_component_files creates files correctly", {
  # Create a temporary directory for output files
  temp_dir <- tempdir()

  # Mock data
  design_name <- "test_design"
  component <- GenomicRanges::GRanges("chr1:100-200")
  names(component) <- "test_component"
  component_name <- "test_components"
  chrom_relative <- GenomicRanges::GRanges("chr1:1-1000")

  # Call the function
  write_component_files(temp_dir, design_name, component, component_name, chrom_relative)

  # Check if files were created
  csv_file <- file.path(temp_dir, "test_design_test_components_1based.csv")
  gff_file <- file.path(temp_dir, "test_design_test_components.gff3")
  expect_true(file.exists(csv_file))
  expect_true(file.exists(gff_file))

  # Clean up the created files
  file.remove(csv_file)
  file.remove(gff_file)
})
