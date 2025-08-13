test_that("get_frame", {
  expect_equal(get_frame(c(10, 5, 12)), c(0, 1, 0))
  expect_equal(get_frame(c(3, 6, 9)), c(0, 0, 0))
  expect_equal(get_frame(c(10, 10, 10)), c(0, 1, 2))
  expect_equal(get_frame(c(11, 11, 11)), c(0, 2, 1))
  expect_equal(get_frame(c(1, 1, 1)), c(0, 1, 2))
  expect_equal(get_frame(c(7, 2, 5, 4)), c(0, 1, 0, 2))
})

test_that("mirror_flip", {
  # 123456789
  #     X
  #  - -
  # flip - -
  ir <- IRanges(start = c(2, 4), width = 1)
  expect_equal(mirror_flip(ir, 4), IRanges(c(8, 6), width = 1))
  # 123456789
  #     X
  #  --
  # flip  --
  ir <- IRanges(start = 2, width = 2)
  expect_equal(mirror_flip(ir, 4), IRanges(7, width = 2))
})

test_that("comb_along", {
  expect_equal(
    sort(comb_along("NN", m = 1, letters = c("A", "B"))),
    c("AN", "BN", "NA", "NB"))
  expect_equal(
    sort(comb_along("NN", m = 2, letters = c("A", "B"))),
    c("AA", "AB", "BA", "BB"))
  expect_equal(
    sort(comb_along("NN", m = 2, letters = c("A"))),
    c("AA"))
  expect_error(
    sort(comb_along("NN", m = 3, letters = c("A"))))
})

test_that("get_genomic_mutation", {
  cds_test <- GRangesList(list("plus" = GRanges(seqnames = "chr",
                                                ranges = IRanges(100, 200),
                                                strand = "+")))
  expect_equal(
    start(get_genomic_mutation(cds_test, 1)),
    start(cds_test[[1]]))
  expect_equal(
    start(get_genomic_mutation(cds_test, 101)),
    end(cds_test[[1]]))

  cds_test <- GRangesList(list("minus" = GRanges(seqnames = "chr",
                                                 ranges = IRanges(c(300, 100), c(400, 200)),
                                                 strand = "-")))
  expect_equal(
    start(get_genomic_mutation(cds_test, 1)),
    end(cds_test[[1]][1]))
  expect_equal(
    start(get_genomic_mutation(cds_test, 202)),
    start(cds_test[[1]][2]))
})
