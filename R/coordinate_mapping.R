#' Create a set of variants that can reverse a mutation injection process.
#'
#' @param variants A GRanges object of the original variants. Must contain
#'   'REF' and 'ALT' metadata columns. Object has to have names.
#' @return A GRanges object representing the reverse variants. The ranges are in
#'   the coordinate system of the *mutated* sequence. The 'ALT' metadata column
#'   contains the sequence needed to revert the mutation (i.e., the original 'REF').
#' @import GenomicRanges
#' @import IRanges
#' @importFrom S4Vectors mcols
reverse_variants <- function(variants) {
  if (is.unsorted(start(variants)) & is.null(names(variants))) {
    stop("Variants must have names or be sorted by start.")
  }
  if (is.unsorted(start(variants))) variants <- variants[order(start(variants))]
  shifts <- nchar(variants$ALT) - nchar(variants$REF)

  if (length(shifts) > 1) {
    # We need the cumulative sum of the shifts *before* the current variant.
    # This is equivalent to lagging the cumulative sum by one.
    cumulative_shifts <- c(0, cumsum(shifts[-length(shifts)]))
  } else if (length(shifts) == 1) {
    cumulative_shifts <- 0
  } else {
    return(GRanges())
  }
  # The region to be replaced in the mutated sequence has a width equal to the
  # length of the ALT allele from the *original* forward variant.
  new_ranges <- IRanges(start = start(variants) + cumulative_shifts,
                        width = nchar(variants$ALT))
  variants_reversed <- GRanges(
    seqnames = seqnames(variants),
    ranges = new_ranges,
    strand = strand(variants),
    REF = variants$ALT,
    ALT = variants$REF
  )
  names(variants_reversed) <- names(variants)
  return(variants_reversed)
}


#' Inject multiple variants and create a coordinate map (Vectorized)
#'
#' This function takes a genomic sequence and a set of variants, injects the
#' ALT alleles, and produces a map to trace any position in the new mutated
#' sequence back to its origin (either the original genome or a specific variant).
#' This version is vectorized for performance and readability.
#'
#' @param variants_in_window A GRanges object of variants.
#' Must have names. Must be sorted.
#' @param seq_len A length of the sequence variants are mapped to.
#' @return A GRanges object mapping the mutated coordinates.
#' @import GenomicRanges
#' @import IRanges
#' @keywords internal
build_variant_layout <- function(variants_in_window, seq_len) {
  if (length(variants_in_window) == 0) {
    return(GRanges(seqnames = "target_seq",
                   ranges = IRanges(1, width = seq_len),
                   source = "genomic",
                   origin_id = "genomic",
                   origin_start = 1,
                   origin_end = seq_len))
  }

  seqlengths(variants_in_window) <- seq_len
  variants_in_window <- sort(variants_in_window)
  if (is.null(names(variants_in_window))) {
    names(variants_in_window) <- seq_along(variants_in_window)
  }

  gap_ranges <- gaps(ranges(variants_in_window), start = 1, end = seq_len)
  gaps_gr <- GRanges(
    seqnames = seqnames(variants_in_window)[1],
    ranges = gap_ranges,
    strand = "*")
  mcols(gaps_gr) <- S4Vectors::DataFrame(REF = NA, ALT = NA, type = "genomic")
  mcols(variants_in_window)$type <- "variant"
  strand(variants_in_window) <- "*"
  all_pieces <- sort(c(gaps_gr, variants_in_window))

  is_variant <- all_pieces$type == "variant"
  mutated_lengths <- width(all_pieces)
  mutated_lengths[is_variant] <- nchar(all_pieces$ALT[is_variant])

  # Calculate the end and start positions in the mutated sequence using cumsum
  mutated_ends <- cumsum(mutated_lengths)
  mutated_starts <- mutated_ends - mutated_lengths + 1

  coordinate_map <- GRanges(
    seqnames = "target_seq",
    ranges = IRanges(start = mutated_starts, end = mutated_ends),
    source = ifelse(is_variant, "variant", "genomic"),
    origin_id = ifelse(is_variant, names(all_pieces), "genomic"),
    origin_start = ifelse(is_variant, 1, start(all_pieces)),
    origin_end = ifelse(is_variant, nchar(all_pieces$ALT), end(all_pieces)))
  return(coordinate_map)
}


#' Describe target composition in genomic and variant coordinates
#'
#' Using a coordinate map, this function translates the location of each target
#' on the mutated sequence into a detailed string describing its constituent
#' genomic and variant parts.
#'
#' @param target A GRanges object with ranges relative to the mutated sequence.
#' @param coordinate_map The map created by `map_variants`.
#' @param window_genomic The GRanges object for the entire genomic window.
#' @return A data.frame containing the original information for each target that
#'   could be successfully mapped, with an added 'coords' metadata column.
#'   The 'coords' column is a semicolon-separated string detailing the genomic
#'   ranges (e.g., 'chr1:100-120') and variant parts (e.g., 'var_rs123:1-10')
#'   that the target covers.
#' @import GenomicRanges
#' @import IRanges
#' @keywords internal
#'
remap_target_to_genomic <- function(
    target, coordinate_map, window_genomic, variants_genomic) {
  if (!methods::is(target, "GRanges")) target <- GRanges(target)

  # If the input is empty, return an empty data.frame with the expected structure
  if (length(target) == 0) {
    return(data.frame(
      seqnames = character(),
      start = integer(),
      end = integer(),
      width = integer(),
      strand = character(),
      coords = character()
    ))
  }

  # Find all overlaps between targets and the coordinate map
  hits <- findOverlaps(target, coordinate_map)

  # Group map hits by the target they overlap
  map_hits_by_target <- split(S4Vectors::subjectHits(hits),
                              S4Vectors::queryHits(hits))

  # Process each target to generate its coordinate string
  results_list <- lapply(seq_along(target), function(i) {
    hit_indices <- map_hits_by_target[[as.character(i)]]

    # If the target does not overlap any part of the map, it cannot be processed.
    if (is.null(hit_indices)) {
      return(NULL) # We'll filter out NULLs later
    }

    target_range <- ranges(target[i])
    map_segments <- coordinate_map[hit_indices]
    # Ensure segments are processed in their order on the mutated sequence
    map_segments <- map_segments[order(start(map_segments))]

    # --- Calculate coordinate components ---
    coord_components <- lapply(seq_along(map_segments), function(j) {
      segment <- map_segments[j]
      overlap_in_mutated <- pintersect(target_range, ranges(segment))

      if (segment$source == "genomic") {
        offset <- start(overlap_in_mutated) - start(segment)
        start_in_window <- segment$origin_start + offset
        end_in_window <- start_in_window + width(overlap_in_mutated) - 1
        # Create a temporary GRanges to easily format the string
        gr_piece <- GRanges(
          seqnames(window_genomic),
          IRanges(
            start = start(window_genomic) + start_in_window - 1,
            end = start(window_genomic) + end_in_window - 1)
        )
        return(as.character(gr_piece))
      } else { # source == "variant"
        offset <- start(overlap_in_mutated) - start(segment)
        start_in_alt <- segment$origin_start + offset
        end_in_alt <- start_in_alt + width(overlap_in_mutated) - 1
        return(paste0(segment$origin_id, ":", start_in_alt, "-", end_in_alt))
      }
    })
    final_coords_str <- paste(coord_components, collapse = ";")
    return(list(coords = final_coords_str))
  })

  # Identify which targets were successfully mapped
  successful_indices <- which(!sapply(results_list, is.null))

  if (length(successful_indices) == 0) {
    return(data.frame(
      seqnames = character(),
      start = integer(),
      end = integer(),
      width = integer(),
      strand = character(),
      coords = character()
    ))
  }

  # Get the results only for the successfully mapped targets
  mapped_results <- results_list[successful_indices]
  coords_df <- do.call(rbind, lapply(mapped_results, as.data.frame))
  original_target_df <- as.data.frame(target[successful_indices])
  colnames(original_target_df) <- colnames(original_target_df)
  original_target_df$strand <- as.character(original_target_df$strand)
  original_target_df$seqnames <- as.character(original_target_df$seqnames)
  final_df <- cbind(original_target_df, coords_df)
  return(final_df)
}
