test_that("get_genomic_seq returns correct DNAStringSet object", {
  seq_result <- get_genomic_seq(
    chromosome = "chr3",
    position = 41224069,
    bp = 20,
    homology_arm_length = 10,
    output_file = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

  expect_s4_class(seq_result, "DNAStringSet")
  expect_equal(length(seq_result), 3)
  expect_equal(nchar(seq_result), c(40, 10, 10))

  expect_equal(
    names(seq_result),
    c("chr3:41224049-41224088:+", "chr3:41224059-41224068:+", "chr3:41224070-41224079:+"))
})

test_that("get_genomic_seq handles positions near chromosome start (trimming)", {
  seq_result <- suppressWarnings(get_genomic_seq(
    chromosome = "chr1",
    position = 5,
    bp = 10,
    homology_arm_length = 3,
    output_file = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  # 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  #           P
  #           1 2 3 4 5  6  7  8 9  10  - Right arm includes position P
  #   4 3 2 1 - Left arm does not include position P
  expect_equal(nchar(seq_result), c(14, 3, 3))
  expect_equal(seq_result[[1]], Biostrings::DNAString(paste(rep("N", 14), collapse="")))
  expect_equal(seq_result[[2]], Biostrings::DNAString("NNN"))
  expect_equal(seq_result[[3]], Biostrings::DNAString("NNN"))

  # Test trimming when left homology arm goes out of bounds
  seq_result <- suppressWarnings(get_genomic_seq(
    chromosome = "chr1",
    position = 1,
    bp = 10,
    homology_arm_length = 3,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))

  expect_equal(nchar(seq_result), c(10, 0, 3))
  expect_equal(seq_result[[1]], Biostrings::DNAString(paste(rep("N", 10), collapse="")))
  expect_equal(seq_result[[2]], Biostrings::DNAString(""))
  expect_equal(seq_result[[3]], Biostrings::DNAString("NNN"))
})

test_that("get_genomic_seq handles positions near chromosome end (trimming)", {
  seq_result <- suppressWarnings(get_genomic_seq(
    chromosome = "chr1",
    position = seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)["chr1"] - 2,
    bp = 10,
    homology_arm_length = 2,
    output_file = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  # 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  #                              P
  #     10 9 8 7 6  5  4  3 2  1 1  2  3- Full seq includes position P
  expect_equal(nchar(seq_result), c(13, 2, 2))
  expect_equal(seq_result[[1]], Biostrings::DNAString(paste(rep("N", 13), collapse="")))
  expect_equal(seq_result[[2]], Biostrings::DNAString("NN"))
  expect_equal(seq_result[[3]], Biostrings::DNAString("NN"))

  seq_result <- suppressWarnings(get_genomic_seq(
    chromosome = "chr1",
    position = seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)["chr1"],
    bp = 10,
    homology_arm_length = 3,
    output_file = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
  # 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  #                                     P
  #           10 9 8 7 6  5  4  3 2  1  - Full seq includes position P
  #                             3 2  1 - Left arm does not include position P
  expect_equal(nchar(seq_result), c(11, 3, 0))
  expect_equal(seq_result[[1]], Biostrings::DNAString(paste(rep("N", 11), collapse="")))
  expect_equal(seq_result[[2]], Biostrings::DNAString("NNN"))
  expect_equal(seq_result[[3]], Biostrings::DNAString(""))
})

test_that("get_genomic_seq uses seqlevelsStyle from genome", {
  seq_result <- get_genomic_seq(
    chromosome = "chr1",
    position = 500,
    bp = 100,
    homology_arm_length = 50,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

    expect_equal(
      seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
      seqlevelsStyle(GRanges(names(seq_result))))
})

test_that("get_genomic_seq handles chromosome not in genome", {
  expect_error(get_genomic_seq(
    chromosome = "chrNonExistent",
    position = 100,
    bp = 10,
    homology_arm_length = 5,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))
})
