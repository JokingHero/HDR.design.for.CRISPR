#' Selects mutations that are non-overlapping at both the genomic and codon level.
#'
#' @param gr A GRanges object of candidate mutations, assumed to be pre-sorted by priority.
#' @param N The maximum number of mutations to select.
#' @return A GRanges object containing the top N non-overlapping mutations.
#' @keywords internal
#'
select_non_overlapping_mutations <- function(gr, N) {
  if (length(gr) == 0) {
    return(gr[0])
  }

  selected_gr <- gr[0]
  seen_codon_keys <- c()
  for (i in seq_along(gr)) {
    if (length(selected_gr) >= N) {
      break
    }

    candidate_mut <- gr[i]
    if (length(findOverlaps(candidate_mut, selected_gr)) > 0) {
      next
    }

    current_cds_df <- candidate_mut$CDS[[1]]
    current_keys <- c()
    has_codon_data <- !is.null(current_cds_df) && nrow(current_cds_df) > 0

    if (has_codon_data) {
      current_keys <- paste(current_cds_df$tx_id, current_cds_df$codon_num, sep = "_")
      if (any(current_keys %in% seen_codon_keys)) {
        next
      }
    }

    selected_gr <- c(selected_gr, candidate_mut)
    if (has_codon_data) {
      seen_codon_keys <- union(seen_codon_keys, current_keys)
    }
  }

  return(selected_gr)
}

#' Calculate compatibility map based on dbSNP information
#' @keywords internal
#'
calculate_compatibility_map <- function(mutations) {
  if (!is.null(mutations$dbSNP)) {
    compat <- sapply(mutations$dbSNP, function(x) {
      if (is.null(x) || nrow(x) == 0) NA else any(x$is_known_variant)
    })
    mutations$compatibility_map <- rep(3, length(mutations)) # Default: not in dbSNP
    mutations$compatibility_map[is.na(compat)] <- 2          # In dbSNP, but not known variant
    mutations$compatibility_map[which(compat)] <- 1          # Known variant
  } else {
    mutations$compatibility_map <- rep(3, length(mutations)) # All unknown
  }
  mutations
}

#' Create the final summary list output
#' @param selected_muts The final GRanges of selected mutations.
#' @param guides The GRanges of guides (can be single or multiple).
#' @param pams The GRanges of pams (can be single or multiple).
#' @param strategy The selection strategy used.
#' @return A list with the standard output format.
#' @keywords internal
#'
format_output_muts <- function(selected_muts, guides, pams, strategy) {
  if (length(selected_muts) == 0) {
    return(list(
      mutations = selected_muts,
      pam_disrupted_count = 0,
      guide_disrupted_count = 0,
      total_cadd = 0,
      any_overlaps_noncoding = FALSE,
      total_compatibility_score = 0
    ))
  }
  strand(selected_muts) <- "+"
  total_cadd <- if (!is.null(selected_muts$CADD)) sum(selected_muts$CADD) else 0
  any_noncoding <- any(selected_muts$overlaps_noncoding > 0)
  sum_compat_map <- sum(selected_muts$compatibility_map)

  if (strategy == "best_global") {
    # For global strategy, count how many of the *original* guides/pams are hit
    pam_disruption_count <- sum(pams %over% selected_muts)
    guide_disruption_count <- sum(guides %over% selected_muts)
  } else {
    # For single-guide strategies, just sum the boolean column
    pam_disruption_count <- sum(selected_muts$pam_disrupted)
    guide_disruption_count <- sum(selected_muts$guide_disrupted)
  }

  list(
    mutations = selected_muts,
    pam_disrupted_count = pam_disruption_count,
    guide_disrupted_count = guide_disruption_count,
    total_cadd = total_cadd,
    any_overlaps_noncoding = any_noncoding,
    total_compatibility_score = sum_compat_map
  )
}

#' @title Select mutation combinations based on different strategies
#' @description
#' A unified function to find optimal or all combinations of mutations.
#' It can operate in three modes:
#' 1.  `"best_single"`: Greedily finds the best set of non-overlapping mutations
#'     for a *single* guide. Prefers mutations that disrupt the PAM, then the guide.
#' 2.  `"best_global"`: Greedily finds the best set of non-overlapping mutations
#'     that disrupt the *maximum number* of guides/PAMs from a larger set.
#' 3.  `"all_valid"`: Finds all possible non-overlapping combinations of a given size
#'     for a *single* guide, filtered to those that overlap the guide or PAM.
#'
#' @param mutations GRanges object with all possible mutations.
#' @param guides A GRanges object. For `strategy = "best_single"` or `"all_valid"`,
#'   this should be a single guide. For `strategy = "best_global"`, this can be
#'   a GRanges object containing multiple guides.
#' @param pams A GRanges object, corresponding to the `guides`.
#' @param maximum_mutations_per_template The number of mutations to select (N).
#' @param strategy A character string specifying the selection mode. One of
#'   `"best_single"`, `"best_global"`, or `"all_valid"`.
#' @return
#' For `strategy = "best_single"` or `"best_global"`, a single list containing
#' the selected mutations and summary statistics.
#' For `strategy = "all_valid"`, a list of lists, where each inner list represents
#' a valid combination.
#' @keywords internal
#'
find_mutation_combinations <- function(
    mutations,
    guides,
    pams,
    maximum_mutations_per_template,
    strategy) {
  strategy <- match.arg(strategy, choices = c("best_single", "best_global", "all_valid"))
  if (length(mutations) == 0) {
    if (strategy == "all_valid") return(list())
    return(format_output_muts(mutations[0], guides, pams, strategy))
  }
  # we set mutations as * to allow overlap with guides/PAMs
  strand(mutations) <- "*"
  mutations <- calculate_compatibility_map(mutations)
  if (is.null(mutations$CADD)) {
    mutations$CADD <- 0
  }
  mutations$overlaps_noncoding <-
    sapply(mutations$noncoding, function(x) if(is.null(x)) 0 else nrow(x))

  if (strategy == "best_global") {
    # Global strategy: counts are used for scoring
    mutations$pam_disrupted <- countOverlaps(mutations, pams)
    mutations$guide_disrupted <- countOverlaps(mutations, guides)
    mutations$distance_to_guide <- sapply(seq_along(mutations),
                                          function(i) min(distance(pams, mutations[i])))
    # Order by most guides/pams disrupted first
    ordering <- order(
      mutations$overlaps_noncoding,
      -mutations$pam_disrupted,
      -mutations$guide_disrupted,
      mutations$CADD,
      mutations$compatibility_map,
      mutations$distance_to_guide
    )
  } else {
    # Single-guide strategies: booleans are used for scoring
    mutations$pam_disrupted <- mutations %over% pams
    mutations$guide_disrupted <- mutations %over% guides
    mutations$distance_to_guide <- distance(mutations, pams)
    # Order by PAM disruption, then guide disruption
    ordering <- order(
      mutations$overlaps_noncoding,
      !mutations$pam_disrupted,
      !mutations$guide_disrupted,
      mutations$CADD,
      mutations$compatibility_map,
      mutations$distance_to_guide
    )
  }
  mutations <- mutations[ordering]

  if (strategy == "all_valid") {
    # Filter to only mutations that actually hit the guide or PAM
    mutations <- mutations[mutations$guide_disrupted | mutations$pam_disrupted]
    if (length(mutations) < maximum_mutations_per_template) return(list())

    # Generate all combinations of the desired size
    combs_indices <- combn(seq_along(mutations), maximum_mutations_per_template, simplify = FALSE)

    result_list <- lapply(combs_indices, function(indices) {
      candidate_muts <- mutations[indices]
      # A combination is valid only if it has no overlapping codons
      selected <- select_non_overlapping_mutations(candidate_muts, maximum_mutations_per_template)
      if (length(selected) == length(candidate_muts)) {
        format_output_muts(selected, guides, pams, strategy)
      } else {
        NA # Mark invalid combinations
      }
    })
    return(result_list[!is.na(result_list)])

  } else { # Handles "best_single" and "best_global"
    # Greedily select the top N non-overlapping mutations
    selected_muts <- select_non_overlapping_mutations(mutations, maximum_mutations_per_template)

    if (strategy == "best_global") {
      # The global strategy also truncates to the exact number requested
      if(length(selected_muts) > maximum_mutations_per_template){
        selected_muts <- selected_muts[1:maximum_mutations_per_template]
      }
    }
    return(format_output_muts(selected_muts, guides, pams, strategy))
  }
}
