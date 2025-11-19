library(testthat)
library(GenomicRanges)
library(IRanges)

test_that("reverse_variants works correctly", {
  variants_std <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(10, 21, 30), end = c(10, 21, 31)),
    strand = "+",
    REF = c("C", "A", "GC"),
    ALT = c("T", "AGGG", "G"))
  names(variants_std) <- c("var_SNV", "var_INS", "var_DEL")
  result <- reverse_variants(variants_std)

  # --- Manual calculation of expected output ---
  # Shifts from original variants:
  # var_SNV: nchar("T") - nchar("C")   = 1 - 1 = 0
  # var_INS: nchar("AGGG") - nchar("A") = 4 - 1 = +3
  # var_DEL: nchar("G") - nchar("GC")  = 1 - 2 = -1

  # Cumulative shifts *before* each variant:
  # For var_SNV: 0
  # For var_INS: 0 (from SNV)
  # For var_DEL: 0 (from SNV) + 3 (from INS) = 3

  # Expected ranges in the mutated sequence:
  # var_SNV_rev: start = 10 + 0 = 10. width = nchar("T") = 1. -> 10-10
  # var_INS_rev: start = 21 + 0 = 21. width = nchar("AGGG")= 4. -> 21-24
  # var_DEL_rev: start = 30 + 3 = 33. width = nchar("G") = 1. -> 33-33
  expected_gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(10, 21, 33), end = c(10, 24, 33)),
    strand = "+",
    REF = c("T", "AGGG", "G"),
    ALT = c("C", "A", "GC"))
  names(expected_gr) <- c("var_SNV", "var_INS", "var_DEL")
  expect_equal(result, expected_gr)

  # Single insertion (positive shift) is handled correctly
  var_single_ins <- variants_std["var_INS"]
  result <- reverse_variants(var_single_ins)
  # Cumulative shift is 0 as it's the first/only variant.
  # Start = 21 + 0 = 21. Width = nchar("AGGG") = 4.
  expected_gr <- GRanges("chr1:21-24", strand = "+", REF = "AGGG", ALT = "A")
  names(expected_gr) <- "var_INS"
  expect_equal(result, expected_gr)

  # Single deletion (negative shift) is handled correctly
  var_single_del <- GRanges("chr1:50-54", REF = "ATCGG", ALT = "A") # del 4bp
  names(var_single_del) <- "big_del"
  result <- reverse_variants(var_single_del)
  # Cumulative shift is 0.
  # Start = 50 + 0 = 50. Width = nchar("A") = 1.
  expected_gr <- GRanges("chr1:50-50", REF = "A", ALT = "ATCGG")
  names(expected_gr) <- "big_del"
  expect_equal(result, expected_gr)

  # Empty input GRanges returns an empty GRanges
  result <- reverse_variants(GRanges())
  expect_equal(result, GRanges())
  expect_s4_class(result, "GRanges")

  # Unsorted variants are handled correctly if they are named
  # Create an unsorted but named version of the standard variants
  variants_unsorted <- variants_std[c(3, 1, 2)]
  result <- reverse_variants(variants_unsorted)
  # The result should be IDENTICAL to the standard case because the function
  # should re-sort them internally before calculating shifts.
  expected_gr <- reverse_variants(variants_std) # Use the sorted version to get expected
  expect_equal(result, expected_gr)

  # Function throws error for unsorted, unnamed variants
  variants_unnamed_unsorted <- variants_std[c(2, 1)]
  names(variants_unnamed_unsorted) <- NULL
  expect_error(
    reverse_variants(variants_unnamed_unsorted),
    "Variants must have names or be sorted by start.")
  result <- reverse_variants(variants_std)
  expect_equal(names(result), names(variants_std))
  expect_equal(result$REF, variants_std$ALT)
  expect_equal(result$ALT, variants_std$REF)
})


test_that("build_variant_layout works properly", {
  # A 100bp sequence with one SNV, one INS, and one DEL.
  seq_len_std <- 100
  variants_std <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(10, 21, 30), end = c(10, 21, 31)),
    strand = "+",
    REF = c("C", "A", "GC"),
    ALT = c("T", "AGGG", "G"))
  result_map <- build_variant_layout(variants_std, seq_len_std)

  # --- Manually calculate the expected output ---
  # Original sequence segments (gaps):
  # 1. 1-9 (length 9)
  # 2. 11-20 (length 10)
  # 3. 22-29 (length 8)
  # 4. 32-100 (length 69)
  # Variant ALT allele lengths: unnamed variants become numeric
  # 1. ("T"): length 1
  # 2. ("AGGG"): length 4
  # 3. ("G"): length 1

  # Mutated sequence layout and lengths:
  # seg1(9) | snv(1) | seg2(10) | ins(4) | seg3(8) | del(1) | seg4(69)
  # Total length = 9 + 1 + 10 + 4 + 8 + 1 + 69 = 102

  # Expected ranges on the target sequence:
  # [1-9], [10-10], [11-20], [21-24], [25-32], [33-33], [34-102]

  expected_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(
      start = c(1, 10, 11, 21, 25, 33, 34),
      end   = c(9, 10, 20, 24, 32, 33, 102)),
    source = c("genomic", "variant", "genomic", "variant", "genomic", "variant", "genomic"),
    origin_id = c("genomic", "1", "genomic", "2", "genomic", "3", "genomic"),
    origin_start = c(1, 1, 11, 1, 22, 1, 32),
    origin_end =   c(9, 1, 20, 4, 29, 1, 100))
  expect_equal(result_map, expected_map)

  # Edge Case: No variants are provided
  result_map <- build_variant_layout(GRanges(), seq_len = 50)
  expected_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(start = 1, end = 50),
    source = "genomic",
    origin_id = "genomic",
    origin_start = 1,
    origin_end = 50)
  expect_equal(result_map, expected_map)
  expect_equal(length(result_map), 1)

  # Edge Case: Variant at the very start of the sequence
  variants_at_start <- GRanges("chr1:1-1", REF = "A", ALT = "T")
  names(variants_at_start) <- "v_start"
  result_map <- build_variant_layout(variants_at_start, seq_len = 50)
  # Expected: a variant segment, then a genomic segment. No leading genomic gap.
  expected_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(start = c(1, 2), end = c(1, 50)),
    source = c("variant", "genomic"),
    origin_id = c("v_start", "genomic"),
    origin_start = c(1, 2),
    origin_end = c(1, 50))
  expect_equal(result_map, expected_map)

  # Edge Case: Variant at the very end of the sequence
  variants_at_end <- GRanges("chr1:50-50", REF = "A", ALT = "TTT") # 3bp INS
  names(variants_at_end) <- "v_end"
  result_map <- build_variant_layout(variants_at_end, seq_len = 50)
  # Expected: a genomic segment, then a variant segment. No trailing genomic gap.
  # Mutated length = 49 (genomic) + 3 (variant) = 52
  expected_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(start = c(1, 50), end = c(49, 52)),
    source = c("genomic", "variant"),
    origin_id = c("genomic", "v_end"),
    origin_start = c(1, 1),
    origin_end = c(49, 3))
  expect_equal(result_map, expected_map)

  # Case: Adjacent variants are handled correctly
  # Two SNVs right next to each other
  adjacent_vars <- GRanges(
    "chr1",
    IRanges(start = c(10, 11), end = c(10, 11)),
    REF = c("A", "G"), ALT = c("T", "C"))
  names(adjacent_vars) <- c("v1", "v2")
  result_map <- build_variant_layout(adjacent_vars, seq_len = 20)
  # Expected: genomic | variant | variant | genomic (no gap between variants)
  # Mutated length = 9 + 1 + 1 + 9 = 20
  expected_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(
      start = c(1, 10, 11, 12),
      end   = c(9, 10, 11, 20)),
    source = c("genomic", "variant", "variant", "genomic"),
    origin_id = c("genomic", "v1", "v2", "genomic"),
    origin_start = c(1, 1, 1, 12),
    origin_end = c(9, 1, 1, 20))
  expect_equal(result_map, expected_map)
})


test_that("remap_target_to_genomic works correclty", {
  # Probe on genomic-only segments is remapped correctly
  # 1. Genomic and Variant Data
  window_genomic <- GRanges("chr1:101-200")
  variants_genomic <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(110, 121, 130), end = c(110, 121, 131)),
    strand = "+",
    REF = c("C", "A", "GC"),
    ALT = c("T", "AGGG", "G"))
  names(variants_genomic) <- c("var_SNV", "var_INS", "var_DEL")

  # 2. The Coordinate Map (manually created for predictability)
  # This map reflects how the target sequence is built from the genomic/variant parts.
  # Target sequence length = 9 (genomic) + 1 (SNV) + 10 (genomic) + 4 (INS) + 8 (genomic) + 1 (DEL) + 69 (genomic) = 102
  coordinate_map <- GRanges(
    "target_seq",
    IRanges(
      start = c(1,  10, 11, 21, 25, 33, 34),
      end   = c(9,  10, 20, 24, 32, 33, 102)),
    source = c("genomic", "variant", "genomic", "variant", "genomic", "variant", "genomic"),
    origin_id = c("genomic", "var_SNV", "genomic", "var_INS", "genomic", "var_DEL", "genomic"),
    origin_start = c(1, 1, 11, 1, 22, 1, 32),
    origin_end =   c(9, 1, 20, 4, 29, 1, 100))
  # extra test for build_variant_layout
  variants_window <- GenomicFeatures::pmapToTranscripts(variants_genomic, window_genomic)
  mcols(variants_window) <- mcols(variants_genomic)
  expect_equal(coordinate_map, build_variant_layout(variants_window, 100))

  # This probe is fully contained within the first genomic segment.
  probes_on_s <- GRanges("target_seq:3-8")
  result <- remap_target_to_genomic(probes_on_s, coordinate_map, window_genomic, variants_genomic)
  # 101 (window start) + 3 (probe start) - 1
  expected <- data.frame(seqnames = "target_seq", start = 3, end = 8, width = 6, strand = "*", coords = "chr1:103-108")
  expect_equal(result, expected)

  # Probe spanning multiple genomic segments is remapped to the union range
  # This probe jumps over the SNV, touching two distinct genomic segments.
  probes_on_s <- GRanges("target_seq:5-15")
  result <- remap_target_to_genomic(probes_on_s, coordinate_map, window_genomic, variants_genomic)
  # The genomic parts are chr1:105-109 and chr1:111-115. The range() of these is 105-115.
  expected <- as.data.frame(probes_on_s)
  expected$strand <- as.character(expected$strand)
  expected$seqnames <- as.character(expected$seqnames)
  expected$coords <- "chr1:105-109;var_SNV:1-1;chr1:111-115"
  expect_equal(result, expected)

  # Probe on a single variant-only segment is remapped to the variant's location
  # This probe is fully contained within the 3bp insertion (which is 4bp on the target seq).
  probes_on_s <- GRanges("target_seq:22-23")
  result <- remap_target_to_genomic(probes_on_s, coordinate_map, window_genomic, variants_genomic)
  expected <- as.data.frame(probes_on_s)
  expected$strand <- as.character(expected$strand)
  expected$seqnames <- as.character(expected$seqnames)
  expected$coords <- "var_INS:2-3"
  expect_equal(result, expected)

  # Probe spanning two variant segments uses the *first* variant's location
  # This probe covers the end of the INS and the start of the DEL.
  probes_on_s <- GRanges("target_seq:24-33:+")
  result <- remap_target_to_genomic(
    probes_on_s, coordinate_map, window_genomic, variants_genomic)
  expected <- as.data.frame(probes_on_s)
  expected$strand <- as.character(expected$strand)
  expected$seqnames <- as.character(expected$seqnames)
  expected$coords <- "var_INS:4-4;chr1:122-129;var_DEL:1-1"
  expect_equal(result, expected)

  # Edge Case: Probe does not overlap any map segments - we should never be here
  gappy_map <- coordinate_map[c(1, 3)]
  probes_on_s <- GRanges("target_seq:10-10") # Probe is in the gap
  result <- remap_target_to_genomic(probes_on_s, gappy_map, window_genomic, variants_genomic)
  expect_true(nrow(result) == 0)

  # Edge Case: Empty inputs are handled gracefully
  probes_on_s <- GRanges()
  result <- remap_target_to_genomic(probes_on_s, coordinate_map, window_genomic, variants_genomic)
  expect_true(nrow(result) == 0)
  expect_s3_class(result, "data.frame")

  # Empty coordinate map
  probes_on_s <- GRanges("target_seq:1-5")
  result_empty_map <- remap_target_to_genomic(probes_on_s, GRanges(), window_genomic, variants_genomic)
  expect_true(nrow(result_empty_map) == 0)

  # Metadata and names are propagated correctly
  probes_on_s <- GRanges("target_seq:3-8", probe_id = "P1", score = 100)
  names(probes_on_s) <- "my_probe_1"
  result <- remap_target_to_genomic(probes_on_s, coordinate_map, window_genomic, variants_genomic)
  expect_equal(rownames(result), "my_probe_1")
  expect_true("probe_id" %in% names(result))
  expect_equal(result$probe_id, "P1")
  expect_equal(result$score, 100)
  expect_true("coords" %in% names(result))

  # Function handles multiple probes correctly
  # A mix of probes to test batch processing
  probes_on_s <- GRanges(
    "target_seq",
    IRanges(start = c(5, 22), end = c(15, 23)),
    probe_type = c("genomic_spanning", "variant_only"))
  names(probes_on_s) <- c("probe1", "probe2")
  result <- remap_target_to_genomic(
    probes_on_s, coordinate_map, window_genomic, variants_genomic)
  expected <- as.data.frame(probes_on_s)
  expected$strand <- as.character(expected$strand)
  expected$seqnames <- as.character(expected$seqnames)
  expected$coords <- c("chr1:105-109;var_SNV:1-1;chr1:111-115",
                       "var_INS:2-3")
  expect_equal(result, expected)
  expect_equal(nrow(result), 2)
  expect_equal(rownames(result), c("probe1", "probe2"))
})
