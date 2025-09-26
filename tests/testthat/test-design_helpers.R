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

test_that("export_design_results creates all files", {
  # Setup a temporary directory
  temp_dir <- tempdir()
  design_name <- "full_export_test"

  # Create mock data required by the function
  mock_variant <- GRanges("chr1:1000-1000", REF = "A", ALT = "G")
  mock_var_data <- GRanges("chr1:950-1050")
  mcols(mock_var_data) <- list(CDS=GRangesList(), dbSNP=GRangesList(), noncoding=GRangesList())
  mock_guides <- GRanges("chr1:980-1000")
  mock_templates <- GRanges("chr1:950-1050", pam_disrupted = 1, guide_disrupted = 1,
                            cadd = 1, snp_quality = 1, overlaps_noncoding=FALSE)
  mock_genome <- BSgenome::BSgenome(
    organism = "Mock sapiens",
    common_name = "Mock",
    provider = "UCSC",
    provider_version = "mock1",
    release_date = "Jan. 2024",
    release_name = "mock1",
    source_url = "http://www.example.com/",
    seqnames = "chr1",
    circ_seqs = character(0),
    mseqnames = NULL,
    seqs_pkgname = "BSgenome.Hsapiens.UCSC.hg38"
  )
  seqlengths(mock_genome) <- c(chr1=2000)
  mockery::stub(export_design_results, 'getSeq',
                mockery::mock(Biostrings::DNAString(paste(rep("A", 2000), collapse=""))))

  # Mock txdb
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(GRanges()))


  # Call the function
  export_design_results(
    output_dir = temp_dir,
    design_name = design_name,
    variant_genomic = mock_variant,
    var_data = mock_var_data,
    guides = mock_guides,
    repair_template = mock_templates,
    genome = mock_genome,
    txdb = txdb
  )

  # Check for expected files
  expected_files <- paste0(design_name, "_", c("query_1based.csv", "query.gff3",
                                               "guides_1based.csv", "guides.gff3",
                                               "templates_1based.csv", "templates.gff3",
                                               "variants_1based.csv", "variants.gff3",
                                               "probes_1based.csv", "probes.gff3"))
  expected_files <- c(paste0(design_name, ".fa"), expected_files)

  # Note: probes are not created if the 'probes_' argument is empty (the default).
  # This test uses the default, so we expect 10 files instead of 12.
  files_exist <- file.exists(file.path(temp_dir, head(expected_files, -2)))
  expect_true(all(files_exist))

  # Clean up created files
  file.remove(file.path(temp_dir, expected_files[files_exist]))
})