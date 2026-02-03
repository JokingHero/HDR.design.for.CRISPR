library(testthat)
library(GenomicRanges)
library(Biostrings)
library(BSgenome)

mock_genome <- DNAStringSet(list(DNAString("AGCTAGCTAG"),
                                 DNAString("TGCATGCATGCA")))
names(mock_genome) <- c("chr1", "chr2")

test_that("is_vcf_invalid correctly identifies rows to fix (Does NOT stop)", {
  variants <- GRanges(
    "chr1",
    IRanges(start = c(2, 3, 4, 5, 6), width = c(1, 2, 1, 3, 1)),
    REF = c("A", "AG", "A", "GTC", ""),
    ALT = c("G", "A", "AT", "", "A")
  )

  # Analysis:
  # 1. A -> G (SNP): Valid.
  # 2. AG -> A (Del, Anchored 'A' matches): Valid.
  # 3. A -> AT (Ins, Anchored 'A' matches): Valid.
  # 4. GTC -> "" (Empty ALT): Invalid.
  # 5. "" -> A (Empty REF): Invalid.

  expect_equal(is_vcf_invalid(variants), c(FALSE, FALSE, FALSE, TRUE, TRUE))

  # The "Rescue" cases (Block Subs / Unanchored)
  variants_rescue <- GRanges(
    "chr1",
    IRanges(start = c(2, 2), width = c(1, 3)),
    REF = c("G", "GCT"),
    ALT = c("TC", "TC")
  )
  # 1. G -> TC (Indel, G != T): Invalid (Unanchored).
  # 2. GCT -> TC (Indel, G != T): Invalid (Unanchored).
  expect_true(all(is_vcf_invalid(variants_rescue)))
})

test_that("assert_vcf_valid acts as the gatekeeper", {
  # Good variants
  good_vars <- GRanges("chr1", IRanges(1,1), REF="A", ALT="G")
  expect_silent(assert_vcf_valid(good_vars))

  # Bad variants (Unfixed)
  bad_vars <- GRanges("chr1", IRanges(1,1), REF="A", ALT="")
  expect_error(assert_vcf_valid(bad_vars))

  bad_vars_2 <- GRanges("chr1", IRanges(1,1), REF="C", ALT="TA")
  expect_error(assert_vcf_valid(bad_vars_2))

  bad_vars_2 <- GRanges("chr1", IRanges(1,1), REF="AC", ALT="T")
  expect_error(assert_vcf_valid(bad_vars_2))
})

test_that("normalize_variants rescues block substitutions", {
  # Genome: A G C T A G ...
  # Index:  1 2 3 4 5 6 ...

  # Scenario: User inputs "C" -> "A" at pos 3.
  # Genome at 3 is C. Preceding base at 2 is G.
  # Result should be: Pos 2, REF="GC", ALT="GA"

  variants <- GRanges(
    "chr1",
    IRanges(start = 3, width = 1),
    REF = "C",
    ALT = "AT"
  )

  normalized <- normalize_variants(variants, mock_genome)

  # Check Assertions passed
  expect_s4_class(normalized, "GRanges")

  # Check Coordinates (Shifted back by 1)
  expect_equal(start(normalized), 2)
  expect_equal(end(normalized), 3) # width 2

  # Check Sequence (Padded with G)
  expect_equal(as.character(normalized$REF), "GC")
  expect_equal(as.character(normalized$ALT), "GAT")
})

test_that("normalize_variants fails gracefully at chromosome start", {
  # Pos 1 is A. We cannot fetch Pos 0.
  # User inputs "A" -> "T" (SNP) -> OK.
  # User inputs "A" -> "TC" (Unanchored Insert) -> FAIL because cannot pad left.

  variants <- GRanges(
    "chr1",
    IRanges(start = 1, width = 1),
    REF = "A",
    ALT = "TC"
  )

  # is_vcf_invalid returns TRUE.
  # correct_variants sees start=1, skips fix.
  # assert_vcf_valid sees invalid variant, throws error.
  expect_error(normalize_variants(variants, mock_genome), "Unanchored indels detected")
})
