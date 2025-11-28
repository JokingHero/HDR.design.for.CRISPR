#' @title Score the sequence even if potentially failing
#' @param sequence A character string representing a DNA sequence.
#' @param score_fun A scoring function that takes a sequence and returns a
#' list/object with a 'score' element.
#' @return The score as a numeric value, or NA if an error occurs.
#' @keywords internal
#'
safely_score_sequence <- function(sequence, score_fun) {
  tryCatch({
    score_fun(sequence)
  }, error = function(e) {
    # On error, return NA. You could also log the error message e.g., warning(e$message)
    return(NA_real_)
  })
}

#' Helper function to find, score, and format guides for one strand
#' This is hardcoded for Cas9 only
#' @keywords internal
process_strand <- function(
    mutated_seq, pam_pattern, strand_char, scorers, score_efficiency,
    padding_score = 35, offset_body = 17, offset_pam = 6) {
    # For Cas9 (default)
    # 23bp guide - 6bp to cut = 17 and 6bp from PAM/start to cut
  is_reverse <- strand_char == "-"
  pam_locs <- IRanges(
    matchPattern(DNAString(pam_pattern), mutated_seq, fixed = FALSE))
  if (length(pam_locs) == 0) return(GRanges())

  # This window is dependent on the outer function
  if (!is_reverse) {
    trim_left  <- padding_score + offset_body
    trim_right <- padding_score + offset_pam
  } else {
    trim_left  <- padding_score + offset_pam
    trim_right <- padding_score + offset_body
  }
  valid_window <- IRanges(
    start = trim_left + 1,
    end   = nchar(mutated_seq) - trim_right
  )
  # --------------------------------------
  # 5'----CUT---NGG        5'----CUT---NGG
  # + 17bp                            -6bp
  # CCN--CUT---3'             CCN--CUT---3'
  # + 6bp                            -17bp
  pam_locs <- pam_locs[pam_locs %over% valid_window]
  if (length(pam_locs) == 0) return(GRanges())

  # Define the full 23bp protospacer + PAM region
  # For NGG (+ strand), the protospacer is upstream (5') of the PAM.
  # For CCN (- strand), the protospacer is downstream (3' on the + strand).
  protospacer_pam_locs <- if (!is_reverse) {
    resize(pam_locs, width = 23, fix = "end")
  } else {
    resize(pam_locs, width = 23, fix = "start")
  }

  guide_locs_20bp <- flank(pam_locs, width = 20, start = !is_reverse)
  guide_seqs_20bp <- extractAt(mutated_seq, guide_locs_20bp)
  if (is_reverse) {
    guide_seqs_20bp <- reverseComplement(guide_seqs_20bp)
  }
  guides_gr <-  GRanges(
    seqnames = "target_seq",
    ranges = protospacer_pam_locs,
    strand = strand_char)
  strand(guides_gr) <- strand_char
  guides_gr$original <- as.character(guide_seqs_20bp)
  guides_gr$with_pam <- as.character(extractAt(mutated_seq, protospacer_pam_locs))

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
#' @param variants_genomic A GRanges object of representing the variants.
#'        Must have an 'ALT' metadata column.
#' @param design_name A character string for naming the design.
#' @param cut_distance_max An integer.
#' @param genome A BSgenome object for sequence retrieval.
#' @param any_ALT_on_guides whether any ALT variant is on the guides
#' @param score_efficiency A logical flag. If TRUE, calculate efficiency scores.
#' @return A GRanges object containing all found guides, their sequences, scores
#'         (for successfully run algorithms), and a final aggregated rank.
#' @keywords internal
#'
get_guides_and_scores <- function(
    variants_genomic, design_name, cut_distance_max, genome,
    any_ALT_on_guides,
    score_efficiency = TRUE) {
  window_width <- width(range(variants_genomic)) +
    (cut_distance_max + 23 - 6) * 2 + 1 + 35 * 2 # 35bp is guide score length padding
  window_genomic <- resize(range(variants_genomic), width = window_width, fix = "center")
  genomic_seq <- suppressWarnings(getSeq(genome, window_genomic))[[1]]

  variants_in_window <- pmapToTranscripts(variants_genomic, window_genomic)
  mcols(variants_in_window) <- mcols(variants_genomic)
  names(variants_in_window) <- names(variants_genomic)

  mutated_seq <- if (any_ALT_on_guides) {
    replaceAt(genomic_seq,
              at = ranges(variants_in_window),
              value = DNAStringSet(variants_genomic$ALT))
  } else { genomic_seq }
  coordinate_map <- if (any_ALT_on_guides) {
    build_variant_layout(variants_in_window, nchar(genomic_seq))
  } else {
    GRanges("target_seq", IRanges(1, nchar(genomic_seq)),
            source = "genomic", origin_id = "genomic",
            origin_start = 1, origin_end = nchar(genomic_seq))
  }

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

  fwd_guides <- process_strand(mutated_seq, "NGG", "+", scorers, score_efficiency)
  rve_guides <- process_strand(mutated_seq, "CCN", "-", scorers, score_efficiency)
  guides <- c(fwd_guides, rve_guides)
  if (length(guides) == 0) {
    # Return an empty GRanges with correct columns if no guides are found
    empty_gr <- GRanges()
    mcols(empty_gr) <- data.frame(
      original = character(), coords = character(), rank_by_scores = numeric())
    if (score_efficiency) {
      for (scorer_name in names(scorers)) {
        mcols(empty_gr)[[scorer_name]] <- numeric()
      }
    }
    return(empty_gr)
  }

  guides <- remap_target_to_genomic(
    guides, coordinate_map, window_genomic, variants_genomic)

  if (score_efficiency) {
    score_cols <- guides[, names(scorers), drop = FALSE]
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
