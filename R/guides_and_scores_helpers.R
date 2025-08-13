#' Safely apply a scoring function to a single sequence
#' @param sequence A character string representing a DNA sequence.
#' @param score_fun A scoring function that takes a sequence and returns a
#' list/object with a 'score' element.
#' @return The score as a numeric value, or NA if an error occurs.
#' @keywords internal
safely_score_sequence <- function(sequence, score_fun) {
  tryCatch({
    score_fun(sequence)
  }, error = function(e) {
    # On error, return NA. You could also log the error message e.g., warning(e$message)
    return(NA_real_)
  })
}

#' Helper function to find, score, and format guides for one strand
#' @keywords internal
process_strand <- function(
    mutated_seq, window_genomic, pam_pattern, strand_char, scorers, score_efficiency) {
  is_reverse <- strand_char == "-"
  pam_locs <- IRanges(matchPattern(DNAString(pam_pattern), mutated_seq))
  if (length(pam_locs) == 0) return(GRanges())

  # This window is dependent on the outer function (35bp of extra for scorers)
  valid_window <- IRanges(start = 36, end = nchar(mutated_seq) - 36)
  pam_locs <- pam_locs[pam_locs %over% valid_window]
  if (length(pam_locs) == 0) return(GRanges())

  guide_locs_20bp <- flank(pam_locs + 1, width = 20, start = !is_reverse)
  guide_seqs_20bp <- extractAt(mutated_seq, guide_locs_20bp)
  if (is_reverse) {
    guide_seqs_20bp <- reverseComplement(guide_seqs_20bp)
  }
  guides_gr <- pmapFromTranscripts(guide_locs_20bp, window_genomic)
  strand(guides_gr) <- strand_char
  guides_gr$original <- as.character(guide_seqs_20bp)

  if (score_efficiency) {
    for (scorer_name in names(scorers)) {
      scorer_info <- scorers[[scorer_name]]
      seqs_for_scoring <- scorer_info$seq_extractor(mutated_seq, pam_locs, is_reverse)
      if (is_reverse) {
        seqs_for_scoring <- as.character(reverseComplement(DNAStringSet(seqs_for_scoring)))
      }
      scores <- sapply(seqs_for_scoring, safely_score_sequence, score_fun = scorer_info$score_fun)
      mcols(guides_gr)[[scorer_name]] <- as.numeric(scores)
    }
  }
  pam_name <- if (is_reverse) "CCN" else "NGG"
  names(guides_gr) <- paste0(pam_name, "_", seq_along(guides_gr))
  return(guides_gr)
}


#' @title Get and Score Guides (Refactored)
#' @description For a given genomic location, this function introduces a mutation,
#' designs potential CRISPR guides, and scores them for efficiency. It is designed
#' to be robust to failures in individual scoring algorithms.
#'
#' @param variant_genomic A GRanges object of length 1 representing the variant.
#'        Must have an 'ALT' metadata column.
#' @param design_name A character string for naming the design.
#' @param guide_distance An integer, the window size (in bp) on each side of the
#'        variant to search for guides.
#' @param genome A BSgenome object for sequence retrieval.
#' @param score_efficiency A logical flag. If TRUE, calculate efficiency scores.
#' @return A GRanges object containing all found guides, their sequences, scores
#'         (for successfully run algorithms), and a final aggregated rank.
#' @keywords internal
#'
get_guides_and_scores_refactored <- function(
    variant_genomic, design_name, guide_distance, genome,
    score_efficiency = TRUE) {
  window_width <- guide_distance * 2 + 1 + 35 * 2 # 35bp is guide score length padding
  window_genomic <- resize(variant_genomic, width = window_width, fix = "center")
  genomic_seq <- getSeq(genome, window_genomic)[[1]]
  mutated_seq <- replaceAt(genomic_seq,
                           at = IRanges(guide_distance + 1 + 35, width = 1),
                           value = DNAStringSet(variant_genomic$ALT))
  scorers <- list(
    doench_2014 = list(
      # Extracts 30bp sequence ending 3bp after PAM
      seq_extractor = function(seq, pams, is_rev) {
        ranges <- resize(pams, width = 27, fix = if (is_rev) "start" else "end")
        ranges <- resize(ranges, width = 30, fix = if (is_rev) "end" else "start")
        as.character(extractAt(seq, ranges))
      },
      score_fun = function(x) crisprScore::getRuleSet1Scores(x)$score
    ),
    # doench_2016 = list(
    #   seq_extractor = function(seq, pams, is_rev) {
    #     ranges <- resize(pams, width = 27, fix = if (is_rev) "start" else "end")
    #     ranges <- resize(ranges, width = 30, fix = if (is_rev) "end" else "start")
    #     as.character(extractAt(seq, ranges))
    #   },
    #   score_fun = function(x) crisprScore::getAzimuthScores(x)$score
    # ),
    # deweirdt_2022 = list(
    #   seq_extractor = function(seq, pams, is_rev) {
    #     ranges <- resize(pams, width = 27, fix = if (is_rev) "start" else "end")
    #     ranges <- resize(ranges, width = 30, fix = if (is_rev) "end" else "start")
    #     as.character(extractAt(seq, ranges))
    #   },
    #   score_fun = function(x) crisprScore::getRuleSet3Scores(x, tracrRNA = "Hsu2013")$score
    # ),
    # wang_2019 = list(
    #   seq_extractor = function(seq, pams, is_rev) {
    #     ranges <- resize(pams, width = 27, fix = if (is_rev) "start" else "end")
    #     as.character(extractAt(seq, ranges))
    #   },
    #   score_fun = function(x) crisprScore::getDeepHFScores(x, enzyme = "WT", promoter = "U6")$score
    # ),
    # kim_2019 = list(
    #   seq_extractor = function(seq, pams, is_rev) {
    #     ranges <- resize(pams, width = 27, fix = if (is_rev) "start" else "end")
    #     ranges <- resize(ranges, width = 30, fix = if (is_rev) "end" else "start")
    #     as.character(extractAt(seq, ranges))
    #   },
    #   score_fun = function(x) crisprScore::getDeepSpCas9Scores(x)$score
    # ),
    moreno_mateos_2015 = list(
      # Extracts 35bp sequence centered on a 29bp core
      seq_extractor = function(seq, pams, is_rev) {
        ranges <- resize(pams, width = 29, fix = if (is_rev) "start" else "end")
        ranges <- resize(ranges, width = 35, fix = if (is_rev) "end" else "start")
        as.character(extractAt(seq, ranges))
      },
      score_fun = function(x) crisprScore::getCRISPRscanScores(x)$score
    ),
    labuhn_2018 = list(
      # Extracts the 20bp protospacer sequence
      seq_extractor = function(seq, pams, is_rev) {
        ranges <- flank(pams + 1, width = 20, start = !is_rev)
        as.character(extractAt(seq, ranges))
      },
      score_fun = function(x) crisprScore::getCRISPRaterScores(x)$score
    )
  )

  fwd_guides <- process_strand(mutated_seq, window_genomic, "GG", "+", scorers, score_efficiency)
  rve_guides <- process_strand(mutated_seq, window_genomic, "CC", "-", scorers, score_efficiency)
  guides <- c(fwd_guides, rve_guides)
  if (length(guides) == 0) {
    # Return an empty GRanges with correct columns if no guides are found
    empty_gr <- GRanges()
    mcols(empty_gr) <- data.frame(original = character(), rank_by_scores = numeric())
    if (score_efficiency) {
      for (scorer_name in names(scorers)) {
        mcols(empty_gr)[[scorer_name]] <- numeric()
      }
    }
    return(empty_gr)
  }

  if (score_efficiency) {
    score_cols <- mcols(guides)[, names(scorers), drop = FALSE]
    successful_scorers <- names(which(sapply(score_cols, function(x) !all(is.na(x)))))

    if (length(successful_scorers) > 0) {
      rank_df <- as.data.frame(
        lapply(score_cols[, successful_scorers, drop=FALSE],
               function(x) rank(-x, na.last = "keep")))
      geom_mean_rank <- apply(rank_df, 1, function(row) {
        valid_ranks <- row[!is.na(row)]
        if (length(valid_ranks) == 0) return(NA)
        exp(mean(log(valid_ranks)))
      })
      guides$rank_by_scores <- rank(geom_mean_rank, na.last = "keep")
    } else {
      guides$rank_by_scores <- NA_real_
    }
  }

  return(guides)
}
