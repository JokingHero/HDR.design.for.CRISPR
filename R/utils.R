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
    gseq,        # DNAString object
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

get_all_possilbe_mutations <- function(
    positions_to_mutate,
    mutation_loci,
    cds_seq) {
  aa_cds_seq <- Biostrings::translate(cds_seq)
  pp <- positions_to_mutate
  # remove positions of our main mutation codon
  this_pos_mutation_loci <- mutation_loci
  codon_count <- ceiling(this_pos_mutation_loci / 3)
  codon_end <- codon_count * 3
  codon_start <- codon_end - 2
  codon_position <- this_pos_mutation_loci - codon_start + 1
  original_codon <- aa_cds_seq[codon_count]
  original_codon_seq <- cds_seq[codon_start:codon_end]

  # 0 is position of the codon_position, therefore
  # filter out other positions of that codon
  if (codon_position == 3) {
    pp <- pp[!pp %in% -2:0]
  } else if (codon_position == 2) {
    pp <- pp[!pp %in% -1:1]
  } else {
    pp <- pp[!pp %in% 0:2]
  }

  # figure out all possible mutations that are synonymous
  sp <- c()
  mutations <- IRanges()
  for (i in pp) {
    this_pos_mutation_loci <- mutation_loci + i
    codon_count <- ceiling(this_pos_mutation_loci / 3)
    codon_end <- codon_count * 3
    codon_start <- codon_end - 2
    codon_position <- this_pos_mutation_loci - codon_start + 1
    original_codon <- aa_cds_seq[codon_count]
    original_codon_seq <- cds_seq[codon_start:codon_end]

    # we want to change only codon_position and keep same codon
    positions_that_we_keep <- c(1:3)[!c(1:3) %in% codon_position]
    original_codon_seq_ <- original_codon_seq[positions_that_we_keep]
    alternate <- Biostrings::GENETIC_CODE[Biostrings::GENETIC_CODE == as.character(original_codon)]

    # filter out original codon from alternate
    alternate <- alternate[names(alternate) != as.character(original_codon_seq)]

    # filter out those codons that don't have the same extension as original codon
    alt_names <- sapply(names(alternate), function(x) {
      paste0(strsplit(x, "")[[1]][positions_that_we_keep], collapse = "")
    })
    alternate <- alternate[alt_names == as.character(original_codon_seq_)]
    if (length(alternate) == 0) next # original position is crucial to this codon
    # instead of picking randomly we will use all available alternate codons here
    replacement <- sapply(names(alternate), function(x) strsplit(x, "")[[1]][codon_position])
    names(replacement) <- NULL
    mut <- IRanges(rep(this_pos_mutation_loci, length(replacement)),
                   width = 1,
                   original = as.character(original_codon_seq[codon_position]),
                   replacement = replacement,
                   shift = i,
                   codon = codon_count)
    mutations <- c(mutations, mut)
  }
  names(mutations) <- seq_along(mutations)
  mutations
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


# mutations can't also be using the same codon more than 1 time
# mutations can't repeat
# we prefer clean mutations (no overlap from annot)
# we prefer mutations that mutate the pam > guide > the rest
# for SNPs we go with 123 as above and for annotations
get_combinations_of_mutations_for_guide <- function(
    mutations, mutations_per_template, pam, guide) {
  mutations$pam_disrupted <- ranges(mutations) %over% ranges(pam)
  mutations$guide_disrupted <- ranges(mutations) %over% ranges(guide)
  mutations$distance_to_guide <- distance(ranges(mutations), ranges(pam))
  mutations$overlaps_something <- mutations$noncoding != "" | mutations$nonsyn_tx_count > 0
  if (!is.null(mutations$compatible)) {
    mutations$compatibility_map <- rep(3, length(mutations))
    mutations$compatibility_map[is.na(mutations$compatible)] <- 2
    mutations$compatibility_map[which(mutations$compatible)] <- 1
  } else {
    mutations$compatibility_map <- rep(2, length(mutations)) # all unknown
  }
  if (is.null(mutations$cadd)) {
    mutations$cadd <- 0
  }
  # FALSES go in front of TRUES
  ordering <- order(mutations$overlaps_something,
                    !mutations$pam_disrupted,
                    !mutations$guide_disrupted,
                    mutations$cadd,
                    mutations$compatibility_map,
                    mutations$distance_to_guide,
                    decreasing = FALSE)
  mutations <- mutations[ordering]
  mutations <- mutations[!duplicated(mutations$codon)] # only one change per codon
  mutations <- mutations[1:mutations_per_template] # might be NULL

  list(mutations, sum(mutations$pam_disrupted), sum(mutations$guide_disrupted), sum(mutations$cadd),
       any(mutations$overlaps_something), sum(mutations$compatibility_map))
}


# mutations can't also be using the same codon more than 1 time
# mutations can't repeat
# we prefer clean mutations (no overlap from annot)
# we prefer mutations that mutate the pam > guide > the rest
# for SNPs we go with 123 as above and for annotations
# we filter only to mutations overlapping guide
get_all_combinations_of_mutations_for_guide <- function(
    mutations, mpt, pam, guide) {
  mutations$pam_disrupted <- ranges(mutations) %over% ranges(pam)
  mutations$guide_disrupted <- ranges(mutations) %over% ranges(guide)
  mutations$distance_to_guide <- distance(ranges(mutations), ranges(pam))
  mutations$overlaps_something <- mutations$noncoding != "" | mutations$nonsyn_tx_count > 0
  if (!is.null(mutations$compatible)) {
    mutations$compatibility_map <- rep(3, length(mutations))
    mutations$compatibility_map[is.na(mutations$compatible)] <- 2
    mutations$compatibility_map[which(mutations$compatible)] <- 1
  } else {
    mutations$compatibility_map <- rep(2, length(mutations)) # all unknown
  }
  if (is.null(mutations$cadd)) {
    mutations$cadd <- 0
  }
  # FALSES go in front of TRUES
  ordering <- order(mutations$overlaps_something,
                    !mutations$pam_disrupted,
                    !mutations$guide_disrupted,
                    mutations$cadd,
                    mutations$compatibility_map,
                    mutations$distance_to_guide,
                    decreasing = FALSE)
  mutations <- mutations[ordering]
  mutations <- mutations[mutations$guide_disrupted | mutations$pam_disrupted, ]
  # now select all possible combinations of N mutations based on the above
  combs <- combn(seq_along(mutations), mpt, simplify = FALSE)
  combs <- combs[sapply(combs, function(x) { # we filter out these combinations that reuse a codon
    length(unique(mutations[x]$codon)) == length(x)
  })]
  lapply(combs, function(x) {
    smut <- mutations[x]
    list(smut, sum(smut$pam_disrupted), sum(smut$guide_disrupted), sum(smut$cadd),
         any(smut$overlaps_something), sum(smut$compatibility_map))
  })
}

#' @title All combinations of mutations
#' @description this is more global version of above
#' we want to find the SNPs that disable as many of the guides as we can
#' we design for single template for all guides
#' @param mutations all possible mutations
#' @param mutations_per_template How many to select?
#' @param pam PAM of the guide
#' @param guide guide
#' @return A list of muts grouped
#' @importFrom utils combn
#'
get_combinations_of_mutations_for_guides <- function(
    mutations, mutations_per_template, pam, guide) {
  mutations$pam_disrupted <- countOverlaps(ranges(mutations), ranges(pam))
  mutations$guide_disrupted <- countOverlaps(ranges(mutations), ranges(guide))
  mutations$overlaps_something <- mutations$noncoding != "" | mutations$nonsyn_tx_count > 0
  if (!is.null(mutations$compatible)) {
    mutations$compatibility_map <- rep(3, length(mutations))
    mutations$compatibility_map[is.na(mutations$compatible)] <- 2
    mutations$compatibility_map[which(mutations$compatible)] <- 1
  } else {
    mutations$compatibility_map <- rep(2, length(mutations)) # all unknown
  }
  mutations$distance_to_guide <- sapply(
    seq_along(mutations),
    function(x) {
      y <- mutations[x]
      min(distance(ranges(pam), ranges(y)))
  })

  # FALSES go in front of TRUES
  ordering <- order(mutations$overlaps_something,
                    -mutations$pam_disrupted,
                    -mutations$guide_disrupted,
                    mutations$compatibility_map,
                    mutations$distance_to_guide,
                    decreasing = FALSE)
  mutations <- mutations[ordering]
  mutations <- mutations[!duplicated(mutations$codon)] # only one change per codon
  mutations <- mutations[1:mutations_per_template] # might be NULL

  list(mutations, sum(ranges(pam) %over% ranges(mutations)),
       sum(ranges(guide) %over% ranges(mutations)),
       any(mutations$overlaps_something), sum(mutations$compatibility_map))
}


annotate_mutations_with_snps <- function(mutations, mutations_genomic, snps) {
  GenomeInfoDb::seqlevelsStyle(mutations_genomic) <- GenomeInfoDb::seqlevelsStyle(snps)
  sbo <- GRanges(BSgenome::snpsByOverlaps(snps, mutations_genomic))
  hits <- findOverlaps(mutations_genomic, sbo)
  mutations$RefSNP_id <- ""
  mutations$alleles_as_ambig <- ""
  mutations$RefSNP_id[S4Vectors::queryHits(hits)] <- sbo$RefSNP_id[S4Vectors::subjectHits(hits)]
  mutations$alleles_as_ambig[S4Vectors::queryHits(hits)] <- sbo$alleles_as_ambig[S4Vectors::subjectHits(hits)]
  mutations$compatible <- mapply(function(base, iupac){
    stringr::str_detect(Biostrings::IUPAC_CODE_MAP[iupac], base)
  }, mutations$replacement, mutations$alleles_as_ambig)
  mutations
}

over_splice_sites <- function(
    mutations_genomic, txdb, intron_bp, exon_bp) {
  ex <- exons(txdb)
  ex <- c(GRanges(seqnames = seqnames(ex),
                  ranges = IRanges(start = start(ex) - intron_bp,
                                   end = start(ex) - exon_bp - 1),
                  strand = "*"),
          GRanges(seqnames = seqnames(ex),
                  ranges = IRanges(start = end(ex) - exon_bp + 1,
                                   end = end(ex) + intron_bp),
                  strand = "*"))
  return(!mutations_genomic %over% ex)
}


annotate_mutations_with_tx <- function(genome, mutations, mutations_genomic, txdb) {
  mutations$nonsyn_tx_count <- 0
  mutations$nonsyn_tx <- ""
  mutations$overlap_tx_count <- 0

  all_tx_cds <- suppressWarnings(GenomicFeatures::cdsBy(
    txdb, by = "tx", use.names = T))
  for (i in seq_along(mutations)) {
    to_check <- all_tx_cds[all_tx_cds %over% mutations_genomic[i]]
    i_cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, to_check)
    i_aa_cds_seq <- suppressWarnings(Biostrings::translate(i_cds_seq))
    i_loc <- GenomicFeatures::mapToTranscripts(mutations_genomic[i], to_check)
    nonsyn <- c()
    for (j in seq_along(to_check)) {
      ji_cds_seq <- i_cds_seq[[j]]
      if (as.character(ji_cds_seq[start(i_loc[j])]) != mutations[i]$original) {
        stop("Nonysynomous SNPs reference mismatch.")
      }
      ji_cds_seq[start(i_loc[j])] <- DNAString(mutations[i]$replacement)
      if (i_aa_cds_seq[[j]] != suppressWarnings(Biostrings::translate(ji_cds_seq))) {
        nonsyn <- c(nonsyn, names(to_check)[j])
      }
    }
    mutations$nonsyn_tx_count[i] <- length(nonsyn)
    mutations$nonsyn_tx[i] <- paste0(nonsyn, collapse = ";")
    mutations$overlap_tx_count[i] <- length(to_check)
  }
  mutations
}

annotate_mutations_with_noncoding <- function(mutations, mutations_genomic, annotation) {
  igff <- rtracklayer::import(annotation)
  igff <- igff[igff$gene_type != "protein_coding"]
  mutations$noncoding <- ""
  for (i in seq_along(mutations)) {
    i_igff <- igff[igff %over% mutations_genomic[i]]
    if (length(i_igff) > 0) {
      mutations$noncoding[i] <-
        paste0(unique(paste0(i_igff$gene_name, " ", i_igff$gene_type)),
               collapse = "; ")
    }
  }
  mutations
}

annotate_mutations_with_cadd <- function(mutations, mutations_genomic, cadd) {
  mutations$cadd <- GenomicScores::gscores(
    cadd, mutations_genomic, ref = mutations$original, alt = mutations$replacement)$default
  mutations
}

melting_temp <- function(x) TmCalculator::Tm_GC(
  ntseq = as.character(x),
  variant = "Primer3Plus", Na = 50, outlist = FALSE)

gc_fract <- function(x) letterFrequency(x, letters = "CG", as.prob = TRUE)

design_probes <- function(mutation_name, st, sp, s, tmin = 59,
                          tmax = 61, len_min = 20, len_max = 25,
                          origin_mut_start = 400 + 1) {
  probes <- GRanges()
  for (len in len_min:len_max) { # for each length
    for (i in st:(sp-len)) { # for each bases
      si <- s[i:(i+len)]
      tm <- melting_temp(si)
      if (tm < tmin | tm > tmax) next
      probes <- c(probes,
                  GRanges(seqnames = mutation_name,
                          ranges = IRanges(start = i, width = len),
                          strand = "+",
                          original = as.character(si),
                          shift = i - origin_mut_start,
                          length = len,
                          Tm = tm,
                          GC = gc_fract(si)))
    }
  }
  return(probes)
}

select_probes <- function(muts_to_cover, candidates) {
  these_probes <- GRanges()
  for (k in 1:length(muts_to_cover)) {
    o <- findOverlaps(candidates, muts_to_cover)
    if (length(S4Vectors::queryHits(o)) == 0) {
      warning("Could not design probes for all of the mutations for ",
              temp_name, " and mutations ", toString(muts_to_cover))
      break
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

get_guides_and_scores <- function(origin_mutation, mutation_name, guide_distance, genomic_seq,
                                  scores = TRUE) {
  window <- resize(origin_mutation, width = guide_distance * 2, fix = "center")
  mutated_seq <- replaceAt(genomic_seq[[1]],
                           at = ranges(origin_mutation),
                           value = DNAStringSet(origin_mutation$replacement))
  pam_fwd <- IRanges(matchPattern(DNAString("GG"), mutated_seq))
  pam_rve <- IRanges(matchPattern(DNAString("CC"), mutated_seq))
  # restrict to window
  pam_fwd <- pam_fwd[pam_fwd %over% ranges(window)]
  pam_rve <- pam_rve[pam_rve %over% ranges(window)]

  # now we make just guides and score them!
  if (length(pam_fwd) > 0) {
    bp27 <- resize(pam_fwd, width = 27, fix = "end")
    bp30 <- as.character(
      extractAt(mutated_seq,
                resize(bp27, width = 30, fix = "start")))
    bp20 <- flank(pam_fwd + 1, width = 20, start = T)

    # doench_2016 <- sapply(bp30, function(x) getAzimuthScores(x)$score) # problems...
    # deweirdt_2022 <- sapply(bp30,
    #   function(x) getRuleSet3Scores(x, tracrRNA = "Hsu2013")$score) # problems
    # wang_2019 <- sapply(as.character(
    #   extractAt(mutated_seq, bp27)),
    #   function(x) getDeepHFScores(x, enzyme = "WT", promoter = "U6")$score) # not windows

    if (scores) {
      doench_2014 <- sapply(bp30, function(x) crisprScore::getRuleSet1Scores(x)$score)
      kim_2019 <- sapply(bp30, function(x) crisprScore::getDeepSpCas9Scores(x)$score)
      moreno_mateos_2015 <- sapply(as.character(
        extractAt(mutated_seq,
                  resize(resize(pam_fwd, width = 29, fix = "end"), width = 35, fix = "start"))),
        function(x) crisprScore::getCRISPRscanScores(x)$score)
      labuhn_2018 <- sapply(as.character(extractAt(mutated_seq, bp20)),
                            function(x) crisprScore::getCRISPRaterScores(x)$score)
    }
    names(bp20) <- paste0("NGG_", as.character(seq_along(bp20)))
    pam_fwd <- GRanges(seqnames = mutation_name,
                       ranges = bp20,
                       strand = "+",
                       original = as.character(extractAt(mutated_seq, bp20)),
                       shift = end(bp20) - start(origin_mutation),
                       doench_2014 = if (scores) doench_2014 else NA,
                       moreno_mateos_2015 = if (scores) moreno_mateos_2015 else NA,
                       labuhn_2018 = if (scores) labuhn_2018 else NA,
                       kim_2019 = if (scores) kim_2019 else NA)
  }

  if (length(pam_rve) > 0) {
    bp27 <- resize(pam_rve, width = 27, fix = "start")
    bp30 <- as.character(
      reverseComplement(extractAt(mutated_seq,
                                  resize(bp27, width = 30, fix = "end"))))
    bp20 <- flank(pam_rve + 1, width = 20, start = F)

    if (scores) {
      doench_2014 <- sapply(bp30, function(x) crisprScore::getRuleSet1Scores(x)$score)
      kim_2019 <- sapply(bp30, function(x) crisprScore::getDeepSpCas9Scores(x)$score)
      moreno_mateos_2015 <- sapply(as.character(
        reverseComplement(extractAt(mutated_seq,
                                    resize(resize(pam_rve, width = 29, fix = "start"), width = 35, fix = "end")))),
        function(x) crisprScore::getCRISPRscanScores(x)$score)
      labuhn_2018 <- sapply(as.character(reverseComplement(extractAt(mutated_seq, bp20))),
                            function(x) crisprScore::getCRISPRaterScores(x)$score)
    }
    names(bp20) <- paste0("CCN_", as.character(seq_along(bp20)))
    pam_rve <- GRanges(seqnames = mutation_name,
                       ranges = bp20,
                       strand = "-",
                       original = as.character(reverseComplement(extractAt(mutated_seq, bp20))),
                       shift = start(bp20) - start(origin_mutation),
                       doench_2014 = if (scores) doench_2014 else NA,
                       moreno_mateos_2015 = if (scores) moreno_mateos_2015 else NA,
                       labuhn_2018 = if (scores) labuhn_2018 else NA,
                       kim_2019 = if (scores) kim_2019 else NA)
  }
  guides <- if (length(pam_fwd) > 0 & length(pam_rve) > 0) {
    c(pam_fwd, pam_rve)
  } else if (length(pam_fwd) > 0) {
    pam_fwd
  } else {
    pam_rve
  }

  if (scores) {
    rank_score <- data.frame(doench_2014 = rank(-1 * guides$doench_2014),
                             moreno_mateos_2015 = rank(-1 * guides$moreno_mateos_2015),
                             labuhn_2018 = rank(-1 * guides$labuhn_2018),
                             kim_2019 = rank(- 1* guides$kim_2019))
    rank_score$rank <- apply(rank_score, 1, function(x) exp(mean(log(x))))
    guides$rank_by_scores <- rank(rank_score$rank)
  } else {
    guides$doench_2014 <- NULL
    guides$moreno_mateos_2015 <- NULL
    guides$labuhn_2018 <- NULL
    guides$kim_2019 <- NULL
  }
  guides
}
