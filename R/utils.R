# translate_safe <- function(dna_string_set) {
#   af <- alphabetFrequency(dna_string_set, baseOnly = FALSE)[, 5:18]
#   ambig <- if (!is.array(af)) {
#     sum(af) > 0
#   } else {
#     rowSums(af) > 0
#   }
#   trans <- suppressWarnings(Biostrings::translate(dna_string_set[!ambig]))
#   final <- rep("NA", length(dna_string_set))
#   final[!ambig] <- as.character(trans)
#   final
# }

translate_safe <- function(dna_string) {
  af <- alphabetFrequency(dna_string, baseOnly = FALSE)[5:18]
  if (sum(af) > 0) {
    "NA"
  } else {
    suppressWarnings(Biostrings::translate(dna_string))
  }
}

translate_robust <- function(dna_string) {
  suppressWarnings(Biostrings::translate(dna_string, if.fuzzy.codon = "solve"))
}


# locus_name = mutation_name
# gseq = genomic_seq[[1]]
# organism_name = organism(genome)
# output_file = file.path(output_dir, paste0(mutation_name, "_cds.gbk"))
# definition = "."
# accession = "."
# version = "."
# keywords = "."
# authors = "snipsnp.com"
# title = "Direct Submission"
# export_source_app = "snipsnp.com"
# comment_text = ""
# mol_type = "genomic DNA"
# division = "UNA"
# codon_start = 1
write_genbank_custom <- function(
    locus_name,         # Character: e.g., "ADA2-203"
    gseq,               # DNAString object
    cds_on_tx,          # IRanges object: CDS exon coordinates relative to gseq
    aa_cds_seq,         # AAString object: translated CDS sequence
    organism_name,      # Character: e.g., "Homo sapiens"
    output_file,        # Character: path for the output .gb file
    definition = ".",
    accession = ".",
    version = ".",
    keywords = ".",
    authors = "snipsnp.com",
    title = "Direct Submission",
    export_source_app = "snipsnp.com", # Tool used for export
    comment_text = "",  # Optional comment line
    mol_type = "genomic DNA",
    division = "UNA",   # GenBank division (e.g., PLN, PRI, HUM, UNA)
    codon_start = 1
) {
  seq_len <- length(gseq)
  date_str <- toupper(format(Sys.Date(), "%d-%b-%Y"))

  # Open file for writing
  file_conn <- file(output_file, "w")
  on.exit(close(file_conn))

  # --- LOCUS Line ---
  # LOCUS       Name      Length bp    mol   topology division Date
  # Example: LOCUS       ADA2-203               40646 bp    DNA     linear   UNA 26-MAY-2025
  # Locus name part is typically padded. Max 16 chars for locus name itself is common.
  locus_line <- sprintf("LOCUS       %-16s %7d bp    DNA     linear   %3s %s",
                        substr(locus_name, 1, 16), # Ensure locus name is not too long for its slot
                        seq_len,
                        division,
                        date_str)
  writeLines(locus_line, file_conn)

  # --- DEFINITION, ACCESSION, VERSION, KEYWORDS ---
  writeLines(paste("DEFINITION ", definition), file_conn)
  writeLines(paste("ACCESSION  ", accession), file_conn) # Note two spaces for ACCESSION
  writeLines(paste("VERSION    ", version), file_conn)   # Note two spaces for VERSION
  writeLines(paste("KEYWORDS   ", keywords), file_conn)  # Note two spaces for KEYWORDS

  # --- SOURCE and ORGANISM ---
  writeLines(paste("SOURCE     ", organism_name), file_conn) # Note one space for SOURCE
  writeLines(paste("  ORGANISM ", organism_name), file_conn)

  # --- REFERENCE ---
  # This is a minimal reference block
  writeLines(sprintf("REFERENCE   1  (bases 1 to %d)", seq_len), file_conn)
  writeLines(paste("  AUTHORS  ", authors), file_conn)
  writeLines(paste("  TITLE    ", title), file_conn)
  journal_line <- sprintf("  JOURNAL   Exported %s from %s", date_str, export_source_app)
  writeLines(journal_line, file_conn)

  # --- COMMENT ---
  if (nzchar(comment_text)) {
    # Handle potential multi-line comments by writing each line prefixed with COMMENT
    comment_lines <- strsplit(comment_text, "\n")[[1]]
    for (cl in comment_lines) {
      writeLines(paste("COMMENT    ", cl), file_conn) # Ensure adequate spacing for COMMENT
    }
  }

  # --- FEATURES ---
  writeLines("FEATURES             Location/Qualifiers", file_conn)

  # Source feature
  writeLines(sprintf("     source          1..%d", seq_len), file_conn)
  writeLines(sprintf("                     /mol_type=\"%s\"", mol_type), file_conn)
  writeLines(sprintf("                     /organism=\"%s\"", organism_name), file_conn)

  # CDS feature
  if (length(cds_on_tx) > 0) {
    cds_loc_strings <- paste0(start(cds_on_tx), "..", end(cds_on_tx))

    # Constructing the join string with wrapping similar to the example
    # Max line width for feature content is around 79 columns.
    # Initial part: "     CDS             join(" (length 26)
    # Continuation indent: "                     " (length 21)

    cds_location_lines <- c()
    # Start with the first part of the join string
    # Max content length on the first feature line (after "     CDS             ")
    max_first_feature_line_content = 79 - nchar("     CDS             ")
    # Max content length on continuation lines (after "                     ")
    max_cont_feature_line_content = 79 - nchar("                     ")

    current_cds_line_content <- paste0("join(", cds_loc_strings[1])

    if (length(cds_loc_strings) > 1) {
      for (j in 2:length(cds_loc_strings)) {
        next_piece_to_add <- paste0(",", cds_loc_strings[j])
        # If this is the first line being built for CDS location
        current_max_len <- if(length(cds_location_lines) == 0) max_first_feature_line_content else max_cont_feature_line_content

        if (nchar(current_cds_line_content) + nchar(next_piece_to_add) <= current_max_len) {
          current_cds_line_content <- paste0(current_cds_line_content, next_piece_to_add)
        } else {
          # Line is full, add it to our collection (with a comma if it's not a natural end)
          cds_location_lines <- c(cds_location_lines, paste0(current_cds_line_content, ",")) # Add comma for continuation
          current_cds_line_content <- cds_loc_strings[j] # Start new line content with this piece
        }
      }
    }
    # Add the final part of current_cds_line_content and the closing parenthesis
    cds_location_lines <- c(cds_location_lines, paste0(current_cds_line_content, ")"))

    # Write the CDS location lines
    writeLines(paste0("     CDS             ", cds_location_lines[1]), file_conn)
    if (length(cds_location_lines) > 1) {
      for (j in 2:length(cds_location_lines)) {
        writeLines(paste0("                     ", cds_location_lines[j]), file_conn)
      }
    }

    # CDS Qualifiers
    writeLines(sprintf("                     /codon_start=%d", codon_start), file_conn)
    writeLines(sprintf("                     /label=\"%s\"", locus_name), file_conn)

    # Translation
    if (length(aa_cds_seq) > 0) {
      if (aa_cds_seq[length(aa_cds_seq)] == AAString("*")) aa_cds_seq <- aa_cds_seq[1:(length(aa_cds_seq)-1)]
      char_aa_seq <- as.character(aa_cds_seq)
      aa_seq_len <- nchar(char_aa_seq)

      translation_output_lines <- c()

      # Max AA chars on the first translation line: 79 - nchar("                     /translation=\"") - 1 (for end quote)
      first_line_aa_max <- 79 - 21 - nchar("/translation=\"") - 1
      # Max AA chars on continuation lines: 79 - nchar("                     ") -1 (for potential end quote)
      cont_line_aa_max <- 79 - 21 -1

      current_pos <- 1
      # First line of translation
      chunk_len <- min(aa_seq_len, first_line_aa_max)
      line_content <- substring(char_aa_seq, current_pos, current_pos + chunk_len - 1)

      if (aa_seq_len <= first_line_aa_max) { # Fits entirely on the first line
        translation_output_lines <- c(translation_output_lines,
                                      sprintf("                     /translation=\"%s\"", line_content))
      } else {
        translation_output_lines <- c(translation_output_lines,
                                      sprintf("                     /translation=\"%s", line_content))
      }
      current_pos <- current_pos + chunk_len

      # Subsequent lines
      while(current_pos <= aa_seq_len) {
        remaining_len <- aa_seq_len - current_pos + 1
        chunk_len <- min(remaining_len, cont_line_aa_max)
        line_content <- substring(char_aa_seq, current_pos, current_pos + chunk_len - 1)

        # Indent subsequent lines
        translation_output_lines <- c(translation_output_lines,
                                      paste0("                     ", line_content))
        current_pos <- current_pos + chunk_len
      }

      # Add closing quote to the very last line of translation
      if (length(translation_output_lines) > 0) {
        translation_output_lines[length(translation_output_lines)] <-
          paste0(translation_output_lines[length(translation_output_lines)], "\"")
      }

      # Write all translation lines
      for(t_line in translation_output_lines){
        writeLines(t_line, file_conn)
      }
    } else { # No amino acid sequence
      writeLines("                     /pseudo", file_conn) # Or /note="coding sequence but no translation provided"
    }
  }

  # --- ORIGIN ---
  writeLines("ORIGIN", file_conn)
  char_genomic_seq <- tolower(as.character(gseq)) # GenBank sequence is often lowercase

  # Sequence formatting: 60 bases per line, 6 blocks of 10, space separated
  # Line numbers are right-justified (width 9)
  for (i in seq(1, seq_len, by = 60)) {
    line_num_str <- sprintf("%9d", i)
    sub_seq_chunk <- substring(char_genomic_seq, i, min(i + 59, seq_len))

    # Split into blocks of 10
    blocks <- c()
    for (j in seq(1, nchar(sub_seq_chunk), by = 10)) {
      blocks <- c(blocks, substring(sub_seq_chunk, j, min(j + 9, nchar(sub_seq_chunk))))
    }
    seq_line_formatted <- paste(blocks, collapse = " ")
    writeLines(paste0(line_num_str, " ", seq_line_formatted), file_conn)
  }

  # --- End of record ---
  writeLines("//", file_conn)

  return(invisible())
}

get_frame <- function(exon_widths) {
  c(0, cumsum(exon_widths) %% 3)[1:length(exon_widths)]
}

#' Flip IRanges around a central mutation point.
#'
#' @param ranges An IRanges object containing the ranges to be flipped, ranges in `sense`
#' @param extension An integer defining the "half-length" of the genomic_seq.
#' The total length will be `2 * extension + 1`, and the mirror point will be `extension + 1`.
#' @return A new IRanges object with the ranges flipped.
#'
mirror_flip <- function(ranges, extension) {
  original_starts <- start(ranges)
  original_ends <- end(ranges)
  flipped_starts <- 2 * (extension + 1) - original_ends
  flipped_ends   <- 2 * (extension + 1) - original_starts
  flipped_iranges <- IRanges(start = flipped_starts, end = flipped_ends)
  return(flipped_iranges)
}


comb_along <- function(seq, m = 2, letters = c("A", "C", "T", "G")) {
  seq <- as.list(strsplit(seq, "")[[1]])
  indices <- utils::combn(seq_along(seq), m)
  letters <- list(letters)
  seq <- apply(indices, 2, function(x) {
    seq[x] <- letters
    do.call(paste0, expand.grid(seq))
  })
  unique(as.vector(seq))
}

get_cds <- function(txdb, ensemble_transcript_id) {
  cds <- suppressWarnings(GenomicFeatures::cdsBy(txdb, by = "tx", use.names = T))
  cds <- cds[ensemble_transcript_id]
  cds <- cds[!duplicated(names(cds))]
  cds
}

get_genomic_mutation <- function(cds, mutation_loci) {
  mut_genomic <- GenomicFeatures::pmapFromTranscripts(
    IRanges(mutation_loci, width = 1), cds)
  mut_genomic <- mut_genomic[[1]][mut_genomic[[1]]$hit]
  mut_genomic
}


# From `mutations` we need to select n groups of `mutations_per_template` mutations
# that can be unique to each template.
# Each codon can't be reused in this calculation.
get_combinations_of_mutations <- function(mutations, n, mutations_per_template) {
  tc <- n
  available_muts <- seq_along(mutations)
  selected_combs <- matrix(nrow = tc, ncol = mutations_per_template)
  i <- 0
  # we just randomize
  stop_counter <- 100000
  stop_counter_i <- 0
  while (i < tc) {
    i <- i + 1
    stop_counter_i <- stop_counter_i + 1
    available_muts <- seq_along(mutations)
    mut123 <- sample(available_muts, size = mutations_per_template, replace = F)
    mut123 <- sort(mut123)
    # mutations can't also be using the same codon more than 1 time
    # mutations can't repeat
    if ((length(unique(mutations[mut123]$codon)) != mutations_per_template) |
        (any(apply(selected_combs, 1, function(x) all(x == mut123)), na.rm = T))) {
      i <- i - 1 # this randomization failed, try again
    } else {
      selected_combs[i, ] <- mut123
    }

    if (stop_counter_i == stop_counter) {
      stop("Can't produce that many templates with these settings.")
    }
  }
  selected_combs
}

#' Vectorized Melting Temperature Calculation (Primer3Plus variant)
#'
#' This function calculates the melting temperature for a DNAStringSet using the
#' TmCalculator::Tm_GC Primer3Plus formula.
#' Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log(Na+)
#' It is significantly faster than calling original Tm_GC.
#'
#' @param x A DNAStringSet object.
#' @param Na The concentration of Na+ in mM.
#' @return A numeric vector of melting temperatures.
melting_temp_vec <- function(x, Na = 50) {
  # Note: We use log10 and convert Na+ to Molar, which is standard.
  gc_prob <- gc_fract(x)
  probe_len <- width(x)

  # The salt correction part of the formula is 16.6 * log10([Na+]) where [Na+] is in Molar.
  # Na is provided in mM, so we divide by 1000.
  salt_correction <- 16.6 * log10(Na / 1000)

  tm <- 81.5 + 0.41 * (gc_prob * 100) - 600 / probe_len + salt_correction
  return(tm)
}

gc_fract <- function(x) letterFrequency(x, letters = "CG", as.prob = TRUE)[, 1]

#' Design Probes on a Sequence and Map them to Genomic Coordinates
#'
#' This function generates candidate probes on a given DNA sequence `s`, filters
#' them by melting temperature, and then uses a coordinate map to accurately

#' translate their positions back to the original genomic coordinates. This
#' approach correctly handles cases where `s` contains indels relative to the
#' original genomic sequence.
#'
#' It also strictly filters out probes where a Variant/SNV overlaps the
#' first 4bp or last 4bp of the probe sequence to ensure stability.
#'
#' @param s A DNAString object representing the sequence to design probes on.
#' @param genomic_context A GRanges object of length 1 representing the original
#'   genomic region from which `s` was derived.
#' @param coordinate_map A GRanges object created by `map_variants`, which describes
#'   the relationship between coordinates on `s` and `genomic_context`.
#' @param variants_genomic A GRanges object of the variants that were injected to
#'   create `s`. This is used by the mapping function to resolve probes that
#'   fall entirely within a variant.
#' @param tmin The minimum desired melting temperature.
#' @param tmax The maximum desired melting temperature.
#' @param len_min The minimum probe length.
#' @param len_max The maximum probe length.
#' @return A GRanges object of valid probes with genomic coordinates and metadata.
#'
design_probes <- function(s, genomic_context, coordinate_map, variants_genomic,
                          tmin = 59, tmax = 61, len_min = 20, len_max = 25) {
  # 1. Generate all possible probe ranges on the input sequence `s`
  all_ranges <- lapply(len_min:len_max, function(len) {
    starts <- 1:(length(s) - len + 1)
    IRanges(start = starts, width = len)
  })
  probe_ranges <- unlist(IRangesList(all_ranges))

  # 2. Extract sequences and calculate melting temperatures
  probe_views <- Views(s, probe_ranges)
  probes_ss <- methods::as(probe_views, "DNAStringSet")
  tms <- melting_temp_vec(probes_ss)

  # 3. Filter probes based on melting temperature
  keep_indices <- which(tms >= tmin & tms <= tmax)
  if (length(keep_indices) == 0) {
    return(GRanges())
  }

  filtered_probes_ss <- probes_ss[keep_indices]
  filtered_ranges <- probe_ranges[keep_indices]
  filtered_tms <- tms[keep_indices]

  # Create a GRanges of the candidate probes relative to 'target_seq'
  probes_on_s <- GRanges("target_seq", ranges = filtered_ranges)
  variant_regions <- coordinate_map[coordinate_map$source == "variant"]
  hits <- findOverlaps(probes_on_s, variant_regions)

  if (length(hits) > 0) {
    p_idx <- queryHits(hits)
    v_idx <- subjectHits(hits)

    # Get the specific ranges involved in the overlap
    p_ranges_hits <- ranges(probes_on_s)[p_idx]
    v_ranges_hits <- ranges(variant_regions)[v_idx]

    # Calculate the intersection (the part of the variant strictly inside the probe)
    overlaps <- pintersect(p_ranges_hits, v_ranges_hits)

    # Calculate the start/end of the variant overlap *relative* to the probe start (1-based)
    # e.g., if probe is 100-120 and variant overlap is 100-101, rel_start is 1
    rel_start <- start(overlaps) - start(p_ranges_hits) + 1
    rel_end   <- end(overlaps) - start(p_ranges_hits) + 1
    probe_lens <- width(p_ranges_hits)

    # Define "Edge": First 4 bp OR Last 4 bp
    # Note: (probe_lens - 3) gives the start of the last 4 bases
    # (e.g., length 20: 17, 18, 19, 20 are the last 4).
    is_edge_overlap <- (rel_start <= 4) | (rel_end >= (probe_lens - 3))

    # Identify indices of probes that fail this check
    bad_probe_indices <- unique(p_idx[is_edge_overlap])

    if (length(bad_probe_indices) > 0) {
      # Remove them from our filtered lists
      keep_mask <- setdiff(seq_along(filtered_ranges), bad_probe_indices)

      if (length(keep_mask) == 0) {
        return(GRanges())
      }

      filtered_probes_ss <- filtered_probes_ss[keep_mask]
      filtered_ranges <- filtered_ranges[keep_mask]
      filtered_tms <- filtered_tms[keep_mask]
      # Re-create the GRanges for the mapping step below
      probes_on_s <- probes_on_s[keep_mask]
    }
  }

  # 4. Map filtered probe ranges from sequence coordinates to genomic coordinates
  # Create a temporary GRanges object for probes relative to sequence `s`
  final_probes <- remap_target_to_genomic(
    target = probes_on_s,
    coordinate_map = coordinate_map,
    window_genomic = genomic_context,
    variants_genomic = variants_genomic)

  # second round of filtering - we want to keep

  # The remap_target_to_genomic function discards metadata, so we re-attach it.
  # The order is preserved.
  final_probes$ALT <- as.character(filtered_probes_ss)
  final_probes$length <- width(filtered_probes_ss)
  final_probes$Tm <- filtered_tms
  final_probes$GC <- gc_fract(filtered_probes_ss)
  return(final_probes)
}

select_probes <- function(muts_to_cover, candidates, temp_name) {
  if (length(muts_to_cover) == 0 || nrow(candidates) == 0) {
    return(data.frame())
  }

  # Keep track of mutations and candidates that are still in play
  muts_remaining <- muts_to_cover
  candidates_remaining <- candidates
  selected_probes_list <- list()

  # --- Greedy Selection Loop ---
  while (length(muts_remaining) > 0) {
    # Get the names of mutations we still need to cover
    mut_names <- names(muts_remaining)

    # 1. Calculate overlaps using string matching
    # For each candidate probe, count how many remaining mutations it covers.
    # We create a matrix where rows are probes, columns are mutations,
    # and values are TRUE/FALSE for coverage. Then, we sum by row.
    # 'fixed = TRUE' makes grepl faster as it does simple string matching.
    overlap_matrix <- sapply(
      mut_names, grepl, candidates_remaining$coords, fixed = TRUE)

    # If the result is a vector (only one mut left), convert it to a matrix
    if (is.vector(overlap_matrix)) {
      overlap_matrix <- matrix(overlap_matrix, ncol = 1)
    }
    overlap_counts <- rowSums(overlap_matrix)

    # Check if any remaining probe can cover any remaining mutation
    if (max(overlap_counts) == 0) {
      warning("Could not design probes for all mutations for ", temp_name,
              ". Uncovered mutations: ", toString(names(muts_remaining)))
      break # Exit the loop, will return probes found so far
    }

    # 2. Find the best probe
    # Order by the number of mutations covered (desc), then by GC content (desc)
    best_probe_idx <- order(overlap_counts, candidates_remaining$GC, decreasing = TRUE)[1]

    # 3. Add the best probe to our results
    best_probe <- candidates_remaining[best_probe_idx, ]
    selected_probes_list[[length(selected_probes_list) + 1]] <- best_probe

    # 4. Update the set of remaining mutations
    # Find which mutations were covered by the selected probe
    # We can reuse the overlap_matrix row for the best probe
    covered_muts_mask <- overlap_matrix[best_probe_idx, ]
    muts_remaining <- muts_remaining[!covered_muts_mask]
    candidates_remaining <- candidates_remaining[-best_probe_idx, ]

    # Safety break: if candidates run out but mutations remain
    if (nrow(candidates_remaining) == 0 && length(muts_remaining) > 0) {
      warning("Ran out of candidate probes for ", temp_name,
              ". Uncovered mutations: ", toString(names(muts_remaining)))
      break
    }
  }

  # --- Finalize and Return ---
  # Combine the list of selected probes into a single data.frame
  if (length(selected_probes_list) > 0) {
    return(do.call(rbind, selected_probes_list))
  } else {
    # Return an empty data.frame with the same structure as candidates
    return(candidates[0, ])
  }
}
