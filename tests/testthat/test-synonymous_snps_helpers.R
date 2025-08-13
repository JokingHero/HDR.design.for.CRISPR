library(testthat)
library(GenomicRanges)

annot <- system.file(
  "tests", "gencode.v42.annotation_mini.gff3", package = "HDR.design.for.CRISPR")
txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(annot)))

test_that("get_search_positions", {
  up_pos <- c(1, 5, 10)
  down_pos <- c(2, 4, 8)
  # --- Test Case 1: Standard SNV on positive strand ---
  # it correctly calculates positions for a positive strand SNV
  variant <- GRanges("chr1:1000-1000", strand = "+")
  result <- get_search_positions(variant, up_pos, down_pos)
  expected_coords <- sort(c(1000 - up_pos, 1000 + down_pos))
  expected_gr <- GRanges("chr1", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)
  expect_equal(start(result), c(990, 995, 999, 1002, 1004, 1008))

  # --- Test Case 2: Standard SNV on negative strand (should be identical) ---
  variant <- GRanges("chr1:1000-1000", strand = "-")
  result <- get_search_positions(variant, up_pos, down_pos)
  # The result should be identical to the positive strand case
  expected_coords <- sort(c(1000 - up_pos, 1000 + down_pos))
  expected_gr <- GRanges("chr1", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)

  # --- Test Case 3: Deletion on positive strand ---
  # A 3-base deletion
  variant <- GRanges("chr5:2000-2002", strand = "+")
  result <- get_search_positions(variant, up_pos, down_pos)
  expected_coords <- sort(c(start(variant) - up_pos, end(variant) + down_pos))
  expected_gr <- GRanges("chr5", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)
  expect_equal(start(result), c(1990, 1995, 1999, 2004, 2006, 2010))
  # --- Test Case 4: Insertion on negative strand ---
  # An insertion is represented by a zero-width range between two bases
  # In GRanges, this is a 1-base range of the anchor base
  variant <- GRanges("chrX:3000-3000", strand = "-")
  result <- get_search_positions(variant, up_pos, down_pos)
  # Logic is the same as for an SNV
  expected_coords <- sort(c(3000 - up_pos, 3000 + down_pos))
  expected_gr <- GRanges("chrX", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)
  expect_equal(start(result), c(2990, 2995, 2999, 3002, 3004, 3008))

  # --- Test Case 5: Only upstream positions provided ---
  variant <- GRanges("chr1:1000-1000", strand = "+")
  result <- get_search_positions(variant, up_pos, c()) # empty downstream
  expected_coords <- sort(1000 - up_pos)
  expected_gr <- GRanges("chr1", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)
  expect_equal(length(result), 3)

  # --- Test Case 6: Only downstream positions provided ---
  variant <- GRanges("chr1:1000-1000", strand = "+")
  result <- get_search_positions(variant, c(), down_pos) # empty upstream
  expected_coords <- sort(1000 + down_pos)
  expected_gr <- GRanges("chr1", IRanges(start = expected_coords, width = 1), strand = "*")
  expect_equal(result, expected_gr)
  expect_equal(length(result), 3)

  # --- Test Case 7: No positions provided ---
  variant <- GRanges("chr1:1000-1000", strand = "+")
  result <- get_search_positions(variant, c(), c())
  expect_equal(result, GRanges())
  expect_equal(length(result), 0)

  # --- Test Case 8: Input validation ---
  # Test with a multi-range GRanges object
  multi_variant <- GRanges("chr1", IRanges(c(100, 200), c(100, 200)))
  expect_error(
    get_search_positions(multi_variant, up_pos, down_pos),
    "`variant_genomic` must be a GRanges object of length 1."
  )
  # Test with non-GRanges input
  expect_error(
    get_search_positions("chr1:100-100", up_pos, down_pos),
    "`variant_genomic` must be a GRanges object of length 1."
  )
})

test_that("get_all_possible_mutations", {
  dummy_genome <- DNAStringSet(c(chr1 = "AGCTGTCA", chr2 = "TTGCA"))
  names(dummy_genome) <- c("chr1", "chr2")
  positions <- GRanges(c("chr1:3:+", "chr2:3:+")) # REF bases are 'C' and 'G'

  result <- get_all_possible_mutations(positions, dummy_genome)

  expect_s4_class(result, "GRanges")
  expect_equal(length(result), 6) # 2 positions * 3 mutations each
  expect_true(all(c("REF", "ALT") %in% names(mcols(result))))
  expect_equal(result$REF, rep(c("C", "G"), each = 3))
  expect_equal(result$ALT, c("A", "T", "G", "A", "C", "T"))
  expect_equal(granges(result), granges(rep(positions, each = 3)))
})

test_that("annotate_variants_with_cds", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  seqlevelsStyle(txdb) <- seqlevelsStyle(genome)

  # plus-strand gene (CTNNB1)
  # This variant is a missense mutation in the CTNNB1 gene on the '+' strand.
  variant_plus <- GRanges("chr3:41224529:+", REF = "A", ALT = "C")
  names(variant_plus) <- "CTNNB1_missense"
  result <- annotate_variants_with_cds(variant_plus, txdb, genome)[[1]]
  expect_s4_class(result, "DataFrame")
  expect_equal(nrow(result), 30)

  target_tx <- "ENST00000396185.8"
  expect_true(target_tx %in% result$tx_id)
  row <- result[result$tx_id == target_tx, ]
  expect_equal(row$tx_id, target_tx)
  expect_equal(row$ranges, 17)
  expect_equal(row$tx_strand, "+")
  expect_equal(row$codon_ref[[1]], "GAT")
  expect_equal(row$codon_alt[[1]], "GCT")
  expect_equal(row$aa_ref[[1]], "D")
  expect_equal(row$aa_alt[[1]], "A")
  expect_equal(row$codon_num, 6)


  # STAT3 minus strand gene
  variant_minus <- GRanges("chr17:42315749:+", REF = "A", ALT = "G")
  names(variant_minus) <- "STAT3_missense"
  result <- annotate_variants_with_cds(variant_minus, txdb, genome)[[1]]
  expect_s4_class(result, "DataFrame")
  expect_equal(nrow(result), 21)

  target_tx <- "ENST00000264657.10"
  expect_true(target_tx %in% result$tx_id)
  row <- result[result$tx_id == target_tx, ]
  expect_equal(row$tx_id, target_tx)
  expect_equal(row$ranges, 2309)
  expect_equal(row$tx_strand, "-")
  expect_equal(row$codon_ref[[1]], "ATG")
  expect_equal(row$codon_alt[[1]], "ACG")
  expect_equal(row$aa_ref[[1]], "M")
  expect_equal(row$aa_alt[[1]], "T")
  expect_equal(row$codon_num, 770)

  # No overlap
  # This variant is on chrX, which is not in our small GFF3 file.
  variant_no_overlap <- GRanges("chrX:1000000:+", REF = "C", ALT = "G")
  names(variant_no_overlap) <- "no_overlap"

  # Run the annotation
  result <- annotate_variants_with_cds(variant_no_overlap, txdb, genome)[[1]]
  expect_s4_class(result, "DataFrame")
  expect_length(result, 0)
  expect_true(nrow(result) == 0)
  expect_true(ncol(result) == 0)
})

test_that("is_outside_splice_sites", {
  # Create a set of variants to test different scenarios
  all_variants_to_test <- GRanges(
    seqnames = c("chr3", "chr3", "chr3", "chr17", "chr17", "chr17", "chr10"),
    ranges = IRanges(
      start = c(
        41223990, # 1. Inside 5' splice site of plus-strand exon
        41224750, # 2. Inside 3' splice site of plus-strand exon
        41224500, # 3. Outside splice site (deep in plus-strand exon)
        42348390, # 4. Inside 5' splice site of minus-strand exon
        42348540, # 5. Inside 3' splice site of minus-strand exon
        42348450, # 6. Outside splice site (deep in minus-strand exon)
        1000      # 7. Completely outside any defined region
      ), width = 1))
  results <- is_outside_splice_sites(
    all_variants = all_variants_to_test,
    txdb = txdb,
    intron_bp = 10,
    exon_bp = 5)
  expect_equal(results, c(F, F, T, F, F, T, T))

  # Test an edge case: a variant right on the border of a splice site region
  variant_on_border <- GRanges("chr3:41223996:+") # Last base of the 5' site
  expect_false(is_outside_splice_sites(variant_on_border, txdb, 10, 5))
  # Test an edge case: a variant just outside the border
  variant_off_border <- GRanges("chr3:41223997:+") # One base after the 5' site
  expect_true(is_outside_splice_sites(variant_off_border, txdb, 10, 5))
})

test_that("annotate_variants_with_snps", {
  our_variants <- GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr1"),
    ranges = IRanges(start = c(100, 200, 300, 400, 100, 500), width = 1),
    strand = "+",
    REF = c("A", "C", "G", "T", "C", "A"),
    ALT = c("C", "G", "T", "A", "G", "T")
  )
  names(our_variants) <- c("multi_match_var", "single_match_var", "allele_mismatch_var",
                           "no_overlap_var", "chr2_var", "na_snp_var")
  known_snps <- GRanges(
    seqnames = c("chr1", "chr1", "chr1", "chr1", "chr2", "chr1"),
    ranges = IRanges(start = c(100, 100, 200, 300, 100, 500), width = 1),
    strand = "+",
    RefSNP_id = c("rs111", "rs222", "rs333", "rs444", "rs555", "rs666"),
    # M = A/C; Y = C/T; R = A/G; S = C/G
    alleles_as_ambig = c("M", "Y", "R", "S", "G", NA))

  results <- annotate_variants_with_snps(our_variants, known_snps)
  expect_s4_class(results, "DataFrameList")
  expect_equal(length(results), length(our_variants))
  expect_equal(names(results), names(our_variants))

  res1 <- results[["multi_match_var"]] # Variant is chr1:100, ALT="C"
  # It should match rs111 (M=A/C) and rs222 (Y=C/T).
  expect_equal(nrow(res1), 2)
  expect_true(all(c("rs111", "rs222") %in% res1$RefSNP_id))
  expect_true(all(res1$is_known_variant)) # 'C' is in 'M' and 'Y'

  res2 <- results[["single_match_var"]] # Variant is chr1:200, ALT="G"
  # It should match rs333 (R=A/G).
  expect_equal(nrow(res2), 1)
  expect_equal(res2$RefSNP_id, "rs333")
  expect_true(res2$is_known_variant) # 'G' is in 'R'

  res3 <- results[["allele_mismatch_var"]] # Variant is chr1:300, ALT="T"
  # It overlaps rs444 (S=C/G), but the allele "T" is not in S.
  expect_equal(nrow(res3), 1)
  expect_equal(res3$RefSNP_id, "rs444")
  expect_false(res3$is_known_variant)

  res4 <- results[["no_overlap_var"]] # Variant is chr1:400
  expect_equal(nrow(res4), 0)

  res5 <- results[["chr2_var"]] # Variant is chr2:100, ALT="G"
  # It should match rs555 (alleles_as_ambig="G").
  expect_equal(nrow(res5), 1)
  expect_equal(res5$RefSNP_id, "rs555")
  expect_true(res5$is_known_variant)

  res6 <- results[["na_snp_var"]] # Variant is chr1:500
  # It overlaps rs666, which has NA for alleles.
  expect_equal(nrow(res6), 1)
  expect_equal(res6$RefSNP_id, "rs666")
  expect_false(res6$is_known_variant)

  empty_variants <- our_variants[0, ]
  results_empty_var <- annotate_variants_with_snps(empty_variants, known_snps)
  expect_s4_class(results_empty_var, "DataFrameList")
  expect_equal(length(results_empty_var), 0)

  empty_snps <- known_snps[0, ]
  results_empty_snp <- annotate_variants_with_snps(our_variants, empty_snps)
  expect_s4_class(results_empty_snp, "DataFrameList")
  expect_equal(length(results_empty_snp), length(our_variants))
  # Check that every resulting DataFrame is empty
  expect_true(all(S4Vectors::isEmpty(results_empty_snp)))
})
