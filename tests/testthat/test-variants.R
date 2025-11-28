test_that("is_vcf_invalid correctly identifies non-VCF compliant variants", {
  variants <- GRanges(
    "chr1",
    IRanges(
      start = 101:105,
      width = c(1, 2, 1, 3, 1), # Width must match nchar(REF)
    ),
    REF = c("A", "AG", "A", "GTC", ""),
    ALT = c("G", "A", "AT", "", "A")
  )

  # Let's analyze each case against VCF rules:
  # 1. ok_snv:             A -> G.  OK.
  # 2. ok_del:             AG -> A. Indel, starts with A. OK.
  # 3. ok_ins:             A -> AT. Indel, starts with A. OK.
  # 4. bad_del_empty:      GTC -> "". VIOLATION: ALT is empty.
  # 5. bad_empty_ref:      "" -> A.  VIOLATION: REF is empty.
  expect_equal(is_vcf_invalid(variants), c(F, F, F, T, T))
  variants <- GRanges(
    "chr1",
    IRanges(
      start = 101:102,
      width = c(3, 1),
    ),
    REF = c("G", "GTC"),
    ALT = c("TC", "TC")
  )
  # 5. bad_ins_unanchored: G -> TC.  VIOLATION: Indel, but G != T.
  # 6. bad_del_unanchored: GTC -> TC. VIOLATION: Indel, but G != T.
  expect_error(is_vcf_invalid(variants))

  # all ok
  variants <- GRanges(
    "chr1",
    IRanges(
      start = c(101, 201, 301),
      width = c(1, 4, 2),
    ),
    REF = c("C", "CTTA", "CA"),
    ALT = c("T", "C", "CAT")
  )
  expect_true(all(!is_vcf_invalid(variants)))

  # empty GRanges
  variants <- GRanges(REF = "", ALT = "")
  expect_equal(length(!is_vcf_invalid(variants)), 0)
})

test_that("reverse_variants handles variants correctly", {
  # SNP
  variants <- GRanges("1", IRanges(10, 10, names = "v1"), REF = "A", ALT = "G")
  rvariants <- reverse_variants(variants)
  expect_equal(
    rvariants,
    GRanges("1", IRanges(10, 10, names = "v1"), REF = "G", ALT = "A")
  )
  # insertion
  variants <- GRanges("1", IRanges(20, 20, names = "v1"), REF = "C", ALT = "CATG")
  rvariants <- reverse_variants(variants)
  expect_equal(
    rvariants,
    GRanges("1", IRanges(20, 23, names = "v1"), REF = "CATG", ALT = "C")
  )
  # deletion
  variants <- GRanges("1", IRanges(30, 33, names = "v1"), REF = "GTCA", ALT = "G")
  rvariants <- reverse_variants(variants)
  expect_equal(
    rvariants,
    GRanges("1", IRanges(30, 30, names = "v1"), REF = "G", ALT = "GTCA")
  )

  # unsorted mix
  dna_seq <- DNAString("AGCTTAGCTAGCTAGCTTAGCTAGCTAGCTTAGCTAGCTCTAGCTTAGCTAGCT")
  variants <- GRanges(
    "chr1",
    IRanges(
      start = c(10, 30, 20),
      end   = c(10, 30, 22),
    ),
    REF = c("A", "C", "TGA"),
    ALT = c("G", "CTTA", "T")
  )
  expect_error(reverse_variants(variants))
  names(variants) <- c("snv", "ins", "del")
  mut_seq <- replaceAt(dna_seq,
    at = ranges(variants),
    value = DNAStringSet(variants$ALT)
  )
  rvariants <- reverse_variants(variants)
  dna_seq2 <- replaceAt(mut_seq,
    at = ranges(rvariants),
    value = DNAStringSet(rvariants$ALT)
  )
  expect_equal(dna_seq2, dna_seq)

  # empty
  variants <- GRanges(REF = character(), ALT = character())
  rvariants <- reverse_variants(variants)
  expect_s4_class(rvariants, "GRanges")
  expect_equal(length(rvariants), 0)
})
