test_that("safely_score_sequence works correctly", {
  # Mock scoring function that succeeds
  successful_scorer <- function(seq) {
    return(list(score = 0.8))
  }

  # Mock scoring function that fails
  failing_scorer <- function(seq) {
    stop("Scoring failed")
  }

  # Test with a successful scorer
  expect_equal(safely_score_sequence("ATCG", successful_scorer), 0.8)

  # Test with a failing scorer
  expect_true(is.na(safely_score_sequence("ATCG", failing_scorer)))
})

# Minimal test for process_strand to check basic guide finding
test_that("process_strand finds guides correctly", {
  # Mock data
  mutated_seq <- Biostrings::DNAString("AGCTAGCTAGCTAGCTAGCTAGCTAGCTGGAGCTAGCT")
  window_genomic <- GenomicRanges::GRanges("chr1:1-40")
  pam_pattern <- "GG"
  strand_char <- "+"
  scorers <- list() # Not testing scoring in this unit test
  score_efficiency <- FALSE

  # Expected result: one guide found
  guides <- process_strand(mutated_seq, window_genomic, pam_pattern, strand_char, scorers, score_efficiency)
  expect_equal(length(guides), 1)
  expect_equal(as.character(guides$original), "GCTAGCTAGCTAGCTAGCTA")
})

# Mock BSgenome object for testing
# This is a simplified mock object for demonstration purposes
# In a real scenario, you might need a more sophisticated mock
# or use a small, real BSgenome object.
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

# Define a DNA sequence for our mock chromosome
seq <- DNAString(paste(rep("ATGC", 500), collapse=""))
# Inject a PAM site for testing
seq[1000:1001] <- "GG"

# Assign the sequence to the mock genome
# This part is tricky as BSgenome objects are not meant to be modified this way.
# A better approach for real tests might be to use a custom BSgenome package
# or to mock the `getSeq` function itself using a library like `mockery`.
# For this example, we'll proceed with a simplified approach, acknowledging its limitations.
# Let's assume we can mock getSeq behavior. For now, we will rely on the fact that
# our test won't actually call the real getSeq but we will mock it.

test_that("get_guides_and_scores_refactored handles no guides found", {
  # Mocking getSeq to return a sequence without PAM sites
  mocker <- mockery::mock(Biostrings::DNAString(paste(rep("AT", 500), collapse="")))
  mockery::stub(get_guides_and_scores_refactored, 'getSeq', mocker)

  variant_genomic <- GRanges("chr1:100:100", REF = "A", ALT = "G")
  guides <- get_guides_and_scores_refactored(variant_genomic, "test", 20, mock_genome, FALSE, FALSE)
  expect_true(length(guides) == 0)
})

test_that("get_guides_and_scores_refactored finds guides correctly", {
  # Mocking getSeq to return a sequence with a PAM site
  mocker <- mockery::mock(Biostrings::DNAString("AGCTAGCTAGCTAGCTAGCTAGCTAGCTGGAGCTAGCT"))
  mockery::stub(get_guides_and_scores_refactored, 'getSeq', mocker)

  variant_genomic <- GRanges("chr1:20:20", REF = "A", ALT = "G")
  guides <- get_guides_and_scores_refactored(variant_genomic, "test", 10, mock_genome, FALSE, FALSE)
  expect_equal(length(guides), 1)
  expect_equal(as.character(guides$original), "GCTAGCTAGCTAGCTAGCTA")
})