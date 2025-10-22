translate_safe <- function(dna_string_set) {
  af <- alphabetFrequency(dna_string_set, baseOnly = FALSE)[, 5:18]
  ambig <- if (!is.array(af)) {
    sum(af) > 0
  } else {
    rowSums(af) > 0
  }
  trans <- suppressWarnings(Biostrings::translate(dna_string_set[!ambig]))
  final <- rep("NA", length(dna_string_set))
  final[!ambig] <- as.character(trans)
  final
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

design_probes <- function(s, template_range,
                          tmin = 59, tmax = 61, len_min = 20, len_max = 25) {
  # 1. Generate all possible probe ranges in a vectorized way
  all_ranges <- lapply(len_min:len_max, function(len) {
    starts <- 1:(nchar(s) - len + 1)
    IRanges(start = starts, width = len)
  })
  probe_ranges <- unlist(IRangesList(all_ranges))

  # 2. Create views on the DNAString to get the probe sequences
  probe_views <- Views(s, probe_ranges)
  probes_ss <- methods::as(probe_views, "DNAStringSet")

  # 3. Calculate melting temperatures for all probes at once
  tms <- melting_temp_vec(probes_ss)

  # 4. Filter probes based on melting temperature
  keep_indices <- which(tms >= tmin & tms <= tmax)

  if (length(keep_indices) == 0) {
    return(GRanges())
  }

  filtered_probes_ss <- probes_ss[keep_indices]
  filtered_ranges <- probe_ranges[keep_indices]
  filtered_tms <- tms[keep_indices]

  # 5. map to genomic coords
  final_probes <- pmapFromTranscripts(filtered_ranges, template_range)

  # 6. Add metadata to the final GRanges object
  final_probes$ALT <- as.character(filtered_probes_ss)
  final_probes$length <- width(filtered_probes_ss)
  final_probes$Tm <- filtered_tms
  final_probes$GC <- gc_fract(filtered_probes_ss)
  return(final_probes)
}

select_probes <- function(muts_to_cover, candidates, temp_name) {
  these_probes <- GRanges()
  strand(muts_to_cover) <- "*"

  for (k in seq_along(muts_to_cover)) {
    o <- findOverlaps(candidates, muts_to_cover)
    if (length(S4Vectors::queryHits(o)) == 0) {
      warning("Could not design probes for all of the mutations for ",
              temp_name, " and mutations ", toString(muts_to_cover))
      return(GRanges())
    }
    o_count <- table(S4Vectors::queryHits(o))
    o_count_ <- rep(0, length(candidates))
    o_count_[as.numeric(names(o_count))] <- o_count
    best_probe <- order(o_count_, candidates$GC, decreasing = T)[1]
    these_probes <- c(these_probes, candidates[best_probe])
    candidates <- candidates[-best_probe]
    muts_to_cover <- muts_to_cover[- S4Vectors::subjectHits(o)[
      S4Vectors::queryHits(o) == best_probe]]
    if (isEmpty(muts_to_cover)) {
      break
    }
  }
  these_probes
}
