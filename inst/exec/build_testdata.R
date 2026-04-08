library(HDR.design.for.CRISPR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(GenomicScores)

# Get temp directory from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # fallback to a default tempdir if not provided
  tmp_base <- tempdir()
} else {
  tmp_base <- args[1]
}

# Load annotation (assuming this script is run from project root)
annot <- "tests/testthat/testdata/gencode.v42.annotation_mini.gff3"

# Rebuild STAT3
output_dir_stat3 <- file.path(tmp_base, "STAT3")
if (!dir.exists(output_dir_stat3)) {
  dir.create(output_dir_stat3, recursive = TRUE)
}

design_hdr(
  design_name = "STAT3",
  optimization_scheme = "balanced",
  seed = 42,
  chrom = "chr17",
  variant_start = 42346635,
  variant_end = 42346635,
  REF = "G",
  ALT = "A",
  ALT_on_genome = TRUE,
  ALT_on_templates = FALSE,
  output_dir = output_dir_stat3,
  annotation = annot,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  maximum_variants_per_template = 3,
  snps = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
  cadd = GenomicScores::getGScores("cadd.v1.6.hg38")
)

# Rebuild CTNNB1
output_dir_ctnnb1 <- file.path(tmp_base, "CTNNB1")
if (!dir.exists(output_dir_ctnnb1)) {
  dir.create(output_dir_ctnnb1, recursive = TRUE)
}

design_hdr(
  design_name = "CTNNB1",
  chrom = "chr3",
  variant_start = 41224622,
  variant_end = 41224622,
  REF = "C",
  ALT = "T",
  ALT_on_genome = TRUE,
  ALT_on_templates = FALSE,
  output_dir = output_dir_ctnnb1,
  annotation = annot,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  maximum_variants_per_template = 3,
  snps = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
  cadd = GenomicScores::getGScores("cadd.v1.6.hg38")
)
