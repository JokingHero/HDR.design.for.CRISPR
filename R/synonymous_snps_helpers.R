#' Generate All Possible Single Nucleotide Variants (SNVs)
#'
#' @description
#' The function returns a `GRanges` object where each input position is expanded
#' to three rows, one for each possible alternate allele. The strand of all
#' input ranges is standardized to `+` to ensure consistent retrieval of the
#' reference base from the forward strand, mimicking VCF format conventions.
#'
#' @param search_positions A `GRanges` object from the `GenomicRanges` package.
#'   Each range should represent a single genomic coordinate (i.e., width of 1).
#' @param genome A `BSgenome` object or a `DNAStringSet` from which the
#'   reference sequence can be extracted using `getSeq()`.
#' @param alphabet A character vector representing the set of possible
#'   nucleotides. Defaults to `c("A", "C", "T", "G")`. The alternate alleles
#'   are selected from this set.
#' @return A `GRanges` object with `length(search_positions) * 3` rows (assuming
#'   the default 4-base alphabet). It includes metadata columns:
#'   \item{REF}{The reference allele at the given position.}
#'   \item{ALT}{The alternate (mutant) allele.}
#' @keywords internal
#'
get_all_possible_mutations <- function(search_positions, genome,
                                       alphabet = c("A", "C", "T", "G")) {
  # because we are treating our variants as in VCF file, we
  # are on plus strand
  strand(search_positions) <- "+"
  REF <- as.character(getSeq(genome, search_positions))
  all_variants <- rep(search_positions, each = (length(alphabet) - 1))
  all_variants$REF <- rep(REF, each = (length(alphabet) - 1))
  all_variants$ALT <- unlist(lapply(REF, function(r) {
    alphabet[alphabet != r]
  }))
  all_variants
}

#' Pre-compute the Integrity State of Transcripts Affected by Variants
#' @keywords internal
#'
analyze_transcript_integrity <- function(
    variants_genomic_on_ts, cds_by_tx, cds_seqs, txdb) {
  # Find which transcripts are affected by ANY template variant.
  variant_hits <- findOverlaps(
    cds_by_tx, variants_genomic_on_ts, ignore.strand = TRUE)
  if (length(variant_hits) == 0) return(list())

  affected_tx_names <- names(cds_by_tx)[S4Vectors::queryHits(variant_hits)]
  exons_by_tx <- suppressWarnings(exonsBy(txdb, by = "tx", use.names = TRUE))

  analysis_list <- lapply(unique(affected_tx_names), function(tx_id) {
    tx_cds_gr <- cds_by_tx[[tx_id]]

    # Identify variants on this specific transcript
    tx_idx_in_cds_list <- which(names(cds_by_tx) == tx_id)
    current_hits_mask <- S4Vectors::queryHits(variant_hits) == tx_idx_in_cds_list
    variants_on_tx_genomic <- variants_genomic_on_ts[S4Vectors::subjectHits(variant_hits)[current_hits_mask]]
    indels_on_tx_genomic <- variants_on_tx_genomic[
      nchar(variants_on_tx_genomic$REF) != nchar(variants_on_tx_genomic$ALT)]

    # Frameshift case. Gene disabled.
    bp_diff <- sum(nchar(variants_on_tx_genomic$REF)) - sum(nchar(variants_on_tx_genomic$ALT))
    if (bp_diff %% 3 != 0) {
      return(list(state = "Frameshifted",
                  broken_from_pos = 1))
    }

    # --- Check for Splice Site disruptions first ---
    if (length(variants_on_tx_genomic) > 0 && tx_id %in% names(exons_by_tx)) {
      tx_all_exons_gr <- exons_by_tx[[tx_id]]
      # we filter to only those exons that overlap CDS
      tx_all_exons_gr <- tx_all_exons_gr[tx_all_exons_gr %over% tx_cds_gr]
      splice_sites <- c(
        resize(tx_all_exons_gr, 2, "start"),
        resize(tx_all_exons_gr, 2, "end"))
      # we ignore strand because variants are in GENOMIC plus strand coords
      # splice sites can be minus strand, but they would be affected
      ss_hits <- findOverlaps(splice_sites, variants_on_tx_genomic, ignore.strand = T)
      if (length(ss_hits) > 0) {
        return(list(state = "Splice site disrupted",
                    broken_from_pos = 1))
      }
    }

    # Inject ALL variants to get the mutated sequence
    variants_on_tx_cds <- pmapToTranscripts(
      variants_on_tx_genomic, cds_by_tx[tx_idx_in_cds_list], ignore.strand = T)
    mcols(variants_on_tx_cds) <- mcols(variants_on_tx_genomic)
    variants_on_tx_cds <- variants_on_tx_cds[order(start(variants_on_tx_cds))]
    original_seq <- cds_seqs[names(cds_seqs) == tx_id][[1]]
    mutated_seq <- replaceAt(
      original_seq,
      ranges(variants_on_tx_cds),
      DNAStringSet(variants_on_tx_cds$ALT))
    # now mutated_seq has potentially shifted coordinates by indels
    indels_on_tx_cds <- variants_on_tx_cds[
      nchar(variants_on_tx_cds$REF) != nchar(variants_on_tx_cds$ALT)]
    coord_map <- if (length(indels_on_tx_cds) > 0) {
      diffs <- nchar(indels_on_tx_cds$ALT) - nchar(indels_on_tx_cds$REF)
      data.frame(
        original_pos = start(indels_on_tx_cds),
        cumulative_diff = cumsum(diffs))
    } else {
      data.frame(
        original_pos = integer(0),
        cumulative_diff = integer(0))
    }

    # --- Early stop codon introduced? ---
    translated_mut_seq <- translate_robust(mutated_seq)
    stop_codons <- which(strsplit(as.character(translated_mut_seq), "")[[1]] == "*")
    # this accounts for all the shifts introduced by the variants
    original_stop_pos <- ceiling(length(mutated_seq) / 3)
    first_ptc_pos <- NA_integer_
    if (length(stop_codons) > 0) {
      premature_stops <- stop_codons[stop_codons < original_stop_pos]
      if (length(premature_stops) > 0) {
        first_ptc_pos <- min(premature_stops)
        return(list(
          state = "PTC introduced",
          broken_from_pos = (first_ptc_pos - 1) * 3 - 1,
          mutated_seq = mutated_seq,
          coord_map = coord_map
        ))
      }
    }

    return(list(
      state = "AA change",
      broken_from_pos = length(mutated_seq), # not broken, last pos
      mutated_seq = mutated_seq,
      coord_map = coord_map
    ))
  })

  names(analysis_list) <- unique(affected_tx_names)
  return(analysis_list)
}

#' Map coordinates from an original sequence to a mutated sequence
#'
#' @param original_positions A numeric vector of positions on the original sequence.
#' @param coordinate_map A data frame created by the new logic, with columns
#'   'original_pos' and 'cumulative_diff'.
#' @return A numeric vector of corresponding positions on the mutated sequence.
map_original_to_mutated_coords <- function(original_positions, coordinate_map) {
  # If there are no indels, the coordinates don't change.
  if (nrow(coordinate_map) == 0) {
    return(original_positions)
  }

  # For each original_position, find the index of the last indel that occurs *before* it.
  # findInterval returns 0 if the position is before the first indel.
  indices <- findInterval(original_positions, coordinate_map$original_pos)

  # Get the corresponding cumulative shifts.
  # If an index is 0, the shift is 0. Otherwise, get it from the map.
  shifts <- ifelse(indices == 0, 0, coordinate_map$cumulative_diff[indices])

  # The new position is the original position plus the total shift from all preceding indels.
  mutated_positions <- original_positions + shifts

  return(mutated_positions)
}

# Helper function to extract codons safely (from your old function)
extract_codon_safe <- function(transcript_seq, codon_num) {
  codon_start_pos <- (codon_num - 1) * 3 + 1
  remaining_bases <- length(transcript_seq) - codon_start_pos + 1

  # Case normal
  if (remaining_bases >= 3) {
    return(subseq(transcript_seq, start = codon_start_pos, width = 3))
  }
  # Case incomplete codon
  if (remaining_bases > 0) {
    partial_codon <- subseq(transcript_seq, start = codon_start_pos, width = remaining_bases)
    padding <- DNAString(strrep("N", 3 - remaining_bases))
    return(xscat(partial_codon, padding))
  }
  # Case nonsense
  return(DNAString("NNN"))
}

#' Annotate SNV effects on CDS, considering background variants.
#'
#' This function annotates a set of single nucleotide variants (`all_variants`)
#' by determining their effect on the coding sequences of overlapping transcripts.
#' It uniquely accounts for a pre-defined set of background variants
#' (`variants_genomic_on_ts`), which can include indels, that modify the
#' transcript sequences before the SNVs of interest are evaluated.
#'
#' @param all_variants A GRanges object of SNVs to be annotated.
#' @param variants_genomic_on_ts A GRanges object of background variants (SNVs or indels)
#'   that define the "mutated" reference transcriptome.
#' @param txdb A TxDb object for transcript information.
#' @param genome A BSgenome object for sequence extraction.
#' @param benign_if_moot_statuses "Frameshifted", "Splice site disrupted", "DOWNSTREAM_OF_PTC"
#'
#' @return A DataFrameList, where each element corresponds to a variant in
#'   `all_variants` and contains a DataFrame of annotations for each affected transcript.
#' @importFrom S4Vectors DataFrame
#' @importFrom methods as
#' @importFrom Biostrings DNAString DNAStringSet translate reverseComplement replaceAt xscat subseq
#' @importFrom GenomicRanges start strand `strand<-`
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment seqnames
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicFeatures mapToTranscripts
#'
annotate_variants_with_cds <- function(
    all_variants, variants_genomic_on_ts, txdb, genome,
    benign_if_moot_statuses = c(
      "Frameshifted",
      "Splice site disrupted",
      "DOWNSTREAM_OF_PTC"
    )) {
  cds_by_tx <- suppressWarnings(cdsBy(txdb, by = "tx", use.names = TRUE))
  cds_seqs <- extractTranscriptSeqs(genome, cds_by_tx)

  # Analyze the impact of background variants on each transcript.
  # This returns the state (e.g., frameshifted, PTC introduced), the new
  # mutated sequence, and a coordinate map to handle indel shifts.
  tx_integrity_info <- analyze_transcript_integrity(
    variants_genomic_on_ts, cds_by_tx, cds_seqs, txdb)

  # --- Annotation Loop ---
  # For each SNV, find its effect on all overlapping transcripts.
  result <- lapply(seq_along(all_variants), function(i) {
    variant <- all_variants[i]
    strand(variant) <- "*"
    hits <- suppressWarnings(mapToTranscripts(
      variant, cds_by_tx, ignore.strand = FALSE))
    if (length(hits) == 0) {
      return(S4Vectors::DataFrame())
    }
    annotations <- S4Vectors::DataFrame()

    # Iterate over each transcript hit for the current variant
    for (j in seq_along(hits)) {
      hit <- hits[j]
      tx_id <- as.character(seqnames(hit))
      tx_strand <- as.character(strand(hit))
      pos_on_original_cds <- start(hit)
      info <- tx_integrity_info[[tx_id]]

      # --- Initialize annotation variables ---
      status <- "SYNONYMOUS" # default status
      pos_on_final_seq <- NA_integer_
      final_seq <- NULL
      codon_num <- NA_integer_
      ref_codon_dna <- DNAString("NNN")
      frame_in_codon <- NA_integer_
      aa_ref <- NA_character_
      aa_alt <- NA_character_

      if (is.null(info)) {
        # Case 1. Unaffected Transcript.
        # No background variants on this transcript.
        # Use original sequence and coords.
        final_seq <- cds_seqs[[tx_id]]
        pos_on_final_seq <- pos_on_original_cds
      } else if (info$state %in% c("Frameshifted", "Splice site disrupted")) {
        # CASE 2: Compromised Transcript.
        # The transcript is fundamentally broken. The SNV's effect is moot.
        status <- info$state
        # All other variables remain NA/default, no codon can be determined.
      } else if (info$state %in% c("PTC introduced", "AA change")) {
        # CASE 3: Modified Transcript.
        # The transcript has a new sequence due to background variants.
        final_seq <- info$mutated_seq
        # Map the SNV's position to the new, potentially shifted, coordinate system.
        pos_on_final_seq <- map_original_to_mutated_coords(pos_on_original_cds, info$coord_map)

        if (info$state == "PTC introduced") {
          # The background variants introduced a Premature Termination Codon.
          # Check if our SNV is upstream or downstream of it.
          if (pos_on_final_seq >= info$broken_from_pos) {
            status <- "DOWNSTREAM_OF_PTC"
            # SNV is in a region that is not translated anymore; no codon calculation.
          } else {
            status <- "AA_CHANGE_BEFORE_PTC"
            # SNV is before the PTC, so its effect can be calculated.
          }
        } else {
          status <- "SYNONYMOUS"
        }
      }

      # --- If a codon can be determined ---

      # This block executes if the status allows for codon calculation
      # (i.e., not compromised or downstream of a PTC).
      if (!is.na(pos_on_final_seq)) {
        codon_num <- floor((pos_on_final_seq - 1) / 3) + 1
        frame_in_codon <- (pos_on_final_seq - 1) %% 3 + 1
        ref_codon_dna <- extract_codon_safe(final_seq, codon_num)

        # Determine effective ALT allele based on transcript strand
        effective_alt <- variant$ALT
        if (tx_strand == "-") {
          effective_alt <- as.character(reverseComplement(DNAString(variant$ALT)))
        }
        alt_codon_dna <- replaceAt(
          ref_codon_dna, IRanges(frame_in_codon, width = 1), effective_alt)

        # Translate to amino acids only if the change is meaningful
        if (status %in% c(
          "", "SYNONYMOUS", "AA_CHANGE", "AA_CHANGE_BEFORE_PTC")) {
          aa_ref <- translate_safe(ref_codon_dna)
          aa_alt <- translate_safe(alt_codon_dna)
        }
      } else {
        alt_codon_dna <- DNAString("NNN")
      }

      tx_annot <- S4Vectors::DataFrame(
        tx_id = tx_id,
        tx_strand = tx_strand,
        pos_on_cds = pos_on_final_seq,
        codon_num = codon_num,
        codon_ref = as.character(ref_codon_dna),
        codon_alt = as.character(alt_codon_dna),
        aa_ref = aa_ref,
        aa_alt = aa_alt,
        status = status)

      # A "disqualifying hit" is a non-synonymous change on a viable transcript.
      is_non_synonymous <- !isTRUE(tx_annot$aa_ref == tx_annot$aa_alt)
      is_on_viable_tx <- !(tx_annot$status %in% benign_if_moot_statuses)
      if (is_non_synonymous && is_on_viable_tx) {
        return(tx_annot)
      }

      annotations <- rbind(annotations, tx_annot)
    }
    return(annotations)
  })
  return(methods::as(result, "DataFrameList"))
}

#' Check positions that are not over splice sites
#' @keywords internal
#'
is_outside_splice_sites <- function(
    all_variants, txdb, intron_bp, exon_bp) {
  ex <- exons(txdb)
  splice_site_regions <- c(
    GRanges(
      seqnames = seqnames(ex),
      ranges = IRanges(start = start(ex) - intron_bp,
                       end   = start(ex) + exon_bp - 1),
      strand = "*"),
    GRanges(
      seqnames = seqnames(ex),
      ranges = IRanges(start = end(ex) - exon_bp + 1,
                       end   = end(ex) + intron_bp),
      strand = "*"))
  splice_site_regions <- reduce(splice_site_regions)
  return(!all_variants %over% splice_site_regions)
}

#' Check if variants are outside CDS start/stop boundaries.
#' @keywords internal
#'
is_outside_cds_boundaries <- function(all_variants, txdb, buffer_bp = 3) {
  cds_by_tx <- cdsBy(txdb, by = "tx", use.names = F)
  if (length(cds_by_tx) == 0) {
    return(rep(TRUE, length(all_variants)))
  }

  cds_spans <- unlist(range(cds_by_tx))
  forbidden_zones <- reduce(c(
    GRanges(
      seqnames = seqnames(cds_spans),
      ranges = IRanges(start = start(cds_spans) - buffer_bp,
                       end   = start(cds_spans) + buffer_bp - 1),
      strand = "*"),
    GRanges(
      seqnames = seqnames(cds_spans),
      ranges = IRanges(start = end(cds_spans) - buffer_bp + 1,
                       end   = end(cds_spans) + buffer_bp),
      strand = "*")))
  return(!all_variants %over% forbidden_zones)
}


#' Check if our variants are knoww synonymous variants
#' @keywords internal
#'
annotate_variants_with_snps <- function(all_variants, snps, alt_allele_col = "ALT") {
  var <- all_variants
  strand(var) <- "*"
  GenomeInfoDb::seqlevelsStyle(var) <- GenomeInfoDb::seqlevelsStyle(snps)[1]
  sbo <- if (methods::is(snps, "ODLT_SNPlocs")) {
    GRanges(BSgenome::snpsByOverlaps(snps, var))
  } else snps

  hits <- findOverlaps(var, sbo)
  if (length(hits) == 0) {
    return(DataFrameList(replicate(length(all_variants), S4Vectors::DataFrame())))
  }

  q_hits <- S4Vectors::queryHits(hits)
  s_hits <- S4Vectors::subjectHits(hits)
  variant_alt_alleles <- mcols(all_variants)[[alt_allele_col]][q_hits]
  snp_info <- mcols(sbo)[s_hits, ]
  all_hits_df <- S4Vectors::DataFrame(
    RefSNP_id = snp_info$RefSNP_id,
    alleles_as_ambig = snp_info$alleles_as_ambig
  )
  all_hits_df$is_known_variant <- mapply(function(base, iupac) {
    if (is.na(iupac) || !iupac %in% names(Biostrings::IUPAC_CODE_MAP)) {
      return(FALSE)
    }
    possible_bases <- Biostrings::IUPAC_CODE_MAP[[iupac]]
    stringr::str_detect(possible_bases, base)
  }, variant_alt_alleles, all_hits_df$alleles_as_ambig)

  grouping_factor <- factor(q_hits, levels = seq_len(length(all_variants)))
  annotations_dfl <- S4Vectors::split(all_hits_df, grouping_factor)

  names(annotations_dfl) <- names(all_variants)
  return(annotations_dfl)
}

#' Check if our variants overlap something non-coding
#' @keywords internal
#'
annotate_variants_with_noncoding <- function(all_variants, annotation) {
  maybe_in_same_dir <- list.files(
    base::dirname(annotation),
    pattern = paste0(tools::file_path_sans_ext(base::basename(annotation)), ".*"))
  is_annot <- tools::file_ext(maybe_in_same_dir) %in% c("gff3", "gtf")
  if (!any(is_annot)) {
    return(DataFrameList(replicate(length(all_variants), S4Vectors::DataFrame())))
  }

  igff <- rtracklayer::import(file.path(base::dirname(annotation), maybe_in_same_dir[is_annot]))
  igff <- igff[igff$transcript_type != "protein_coding" &
                 igff$gene_type != "protein_coding" &
                 igff$type == "transcript", ]
  noncoding <- DataFrameList()
  strand(all_variants) <- "*"
  for (i in seq_along(all_variants)) {
    i_igff <- igff[igff %over% all_variants[i]]
    noncoding[[i]] <- if (length(i_igff) > 0) {
        dt <- S4Vectors::DataFrame(
          tx_id = i_igff$transcript_id,
          gene_name = i_igff$gene_name,
          gene_type = i_igff$gene_type)
        dt <- dt[!duplicated(dt$tx_id) & !is.na(dt$tx_id), ]
        dt
    } else {
      S4Vectors::DataFrame()
    }
  }
  names(noncoding) <- names(all_variants)
  return(noncoding)
}

#' CADD score
#' @keywords internal
#'
annotate_variants_with_cadd <- function(all_variants, cadd) {
  GenomicScores::gscores(
    cadd, all_variants,
    ref = all_variants$REF,
    alt = all_variants$ALT)$default
}

#' Annotate Variants with AlphaGenome Scores
#'
#' This function takes a GRanges object containing genetic variants, calls a
#' Python script to score them using the AlphaGenome API, and returns the
#' the new scores as DataFrameList.
#' @param all_variants A \code{GRanges} object. Must contain metadata columns
#'   named \code{REF} and \code{ALT} for the reference and alternate alleles.
#' @param alphagenome_key A character string containing your private AlphaGenome API key.
#' @param species A character string, either "human", "mouse".
#' @param python_exec The command to execute Python (e.g., "python", "python3").
#'   Default is "python".
#' @return `all_variants` with alphagenome mcols
#' @importFrom rlang .data
#' @keywords internal
#'
annotate_mutations_with_alphagenome <- function(all_variants,
                                                alphagenome_key,
                                                species,
                                                python_exec = "python3",
                                                alphagenome_context = "") {
  if (Sys.which(python_exec) == "") {
    stop(paste0("Python executable '", python_exec, "' not found. ",
                "Please make sure python is installed and in your PATH."))
  }
  species <- match.arg(species, c("human", "mouse"))
  script_path <- system.file(
    "exec", "score_variants_with_alphagenome.py", package = "HDR.design.for.CRISPR")
  if (script_path == "") {
    stop(paste0(
      "Could not find 'score_variants_with_alphagenome.py' in the 'exec' folder of the package. Please ensure the script is in the 'inst/exec' directory and re-install the package."
    ))
  }

  temp_input_csv <- tempfile(fileext = ".csv")
  temp_output_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(temp_input_csv, temp_output_csv), force = TRUE))
  variants_df <- as.data.frame(all_variants)
  variants_df <- dplyr::select(variants_df, "seqnames", "start", "REF", "ALT")
  readr::write_csv(variants_df, temp_input_csv)
  args <- c(
    shQuote(script_path),
    "--api_key", shQuote(alphagenome_key),
    "--input", shQuote(temp_input_csv),
    "--species", species,
    "--output", shQuote(temp_output_csv)
  )
  result <- system2(python_exec, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop(
      "Python script execution failed with exit status ", status, ".\n",
      "Command: ", python_exec, " ", paste(args, collapse = " "), "\n",
      "Output:\n", paste(result, collapse = "\n")
    )
  }

  if (!file.exists(temp_output_csv) || file.info(temp_output_csv)$size == 0) {
    stop("Python script seems to have completed, but the output file is missing or empty.\n",
             "Output:\n", paste(result, collapse = "\n"))
  }
  quantile_score <- variant_id <- output_type <- max_score <-
    alphagenome_sites_max <- alphagenome_usage_max <-
    alphagenome_junctions_max <- alphagenome_composite_score <- NULL
  ag_dt <- readr::read_csv(temp_output_csv, show_col_types = FALSE)
  required_input_cols <- c("variant_id", "output_type", "quantile_score")
  if (!all(required_input_cols %in% names(ag_dt))) {
    stop(
      "AlphaGenome output is missing required columns: ",
      paste(setdiff(required_input_cols, names(ag_dt)), collapse = ", ")
    )
  }
  # 1. Take Absolute Value (Magnitude of effect)
  # We care about distance from median (0), whether gain or loss of splicing.
  ag_dt$quantile_score <- abs(ag_dt$quantile_score)

  has_context <- !is.null(alphagenome_context) && alphagenome_context != ""
  if (has_context && "biosample_name" %in% names(ag_dt)) {
    # If user wants specific tissue:
    # - KEEP rows matching that tissue.
    # - ALSO KEEP 'SPLICE_SITES'. Sites are based on DNA motifs and often don't have
    #   tissue labels, but a broken motif is broken everywhere.
    ag_dt <- ag_dt[
      (ag_dt$biosample_name == alphagenome_context) |
        (ag_dt$output_type == "SPLICE_SITES") |
        is.na(ag_dt$biosample_name),
    ]
  }
  # If NO context (Global): We keep ALL rows to find the "worst-case" tissue.
  ag_dt <- dplyr::group_by(ag_dt, variant_id, output_type)
  ag_dt <- suppressMessages(dplyr::summarise(
    ag_dt,
    max_score = max(quantile_score, na.rm = TRUE),
    .groups = "drop"
  ))
  ag_dt$max_score[!is.finite(ag_dt$max_score)] <- 0

  ag_wide <- tidyr::pivot_wider(
    ag_dt,
    id_cols = variant_id,
    names_from = output_type,
    values_from = max_score,
    values_fill = 0)

  required_cols <- c("SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
  for (col in required_cols) {
    if (!col %in% names(ag_wide)) ag_wide[[col]] <- 0
  }

  # Formula: Sites + Usage + (Junctions / 5)
  # Junctions are normalized by 5 because raw junction counts have a larger dynamic range.
  ag_wide$alphagenome_composite_score <- (
    ag_wide$SPLICE_SITES +
      ag_wide$SPLICE_SITE_USAGE +
      (ag_wide$SPLICE_JUNCTIONS / 5.0)
  )
  ag_wide$alphagenome_sites_max <- ag_wide$SPLICE_SITES
  ag_wide$alphagenome_usage_max <- ag_wide$SPLICE_SITE_USAGE
  ag_wide$alphagenome_junctions_max <- ag_wide$SPLICE_JUNCTIONS
  ag_wide <- dplyr::select(
    ag_wide,
    variant_id,
    alphagenome_sites_max,
    alphagenome_usage_max,
    alphagenome_junctions_max,
    alphagenome_composite_score
  )

  original_ids_order <- paste0(
    as.character(GenomicRanges::seqnames(all_variants)), ":",
    GenomicRanges::start(all_variants), ":",
    all_variants$REF, ">",
    all_variants$ALT)

  out <- data.frame(variant_id = original_ids_order, stringsAsFactors = FALSE)
  out <- dplyr::left_join(out, ag_wide, by = "variant_id")
  out$variant_id <- NULL
  for (col in names(out)) {
    if (is.numeric(out[[col]])) out[[col]][is.na(out[[col]])] <- 0
  }
  return(out)
}

#' Prepare Candidate Synonymous Mutations
#' @description Generates, filters, and annotates all possible synonymous SNPs.
#' @keywords internal
prepare_candidate_snps <- function(
    candidate_snp_map,
    annotation, txdb, genome,
    variants_genomic_on_ts,
    intron_bp, exon_bp, clinvar, snps, cadd,
    alphagenome_key, python_exec, alphagenome_context) {

  # some variants might be indels - complicates codon assesments!
  is_indel <- nchar(variants_genomic_on_ts$REF) != nchar(variants_genomic_on_ts$ALT)
  search_positions <- GRanges(candidate_snp_map$coords)
  search_positions <- sort(unique(search_positions))
  # remove indel variant positions - but this is already done trough the candidate_snp_map
  search_positions <- search_positions[
    !(search_positions %over% variants_genomic_on_ts[!is_indel])]
  mcols(search_positions) <- NULL
  # special case where there is nothing to design in terms of SNPs
  # all guides are knocked out with its own variant?
  if (length(search_positions) == 0) {
    message("No valid positions found for introducing synonymous mutations.")
    return(GRanges())
  }
  # Generate all possible 3 SNVs for each unique position.
  all_variants <- get_all_possible_mutations(search_positions, genome)

  # filter out too close to the splice sites
  is_outside <- is_outside_splice_sites(all_variants, txdb, intron_bp, exon_bp)
  all_variants <- all_variants[is_outside]
  if (length(all_variants) == 0) return(GRanges())
  # filter out variants near the start/stop of coding sequences
  is_outside_cds <- is_outside_cds_boundaries(all_variants, txdb, 3)
  all_variants <- all_variants[is_outside_cds]
  if (length(all_variants) == 0) return(GRanges())

  # Filter out known pathogenic ClinVar variants
  if (!is.null(clinvar)) {
    # TODO This logic should be improved to check for specific ALT alleles
    # and pathogenicity, but for now, we just exclude any overlapping position.
    vcf <- VariantAnnotation::readVcf(clinvar)
    vcf <- SummarizedExperiment::rowRanges(vcf)
    vcf <- vcf[lengths(vcf) == 1]
    seqlevelsStyle(vcf) <- seqlevelsStyle(all_variants)
    is_not_over_clinvar <- !(all_variants %over% vcf)
    all_variants <- all_variants[is_not_over_clinvar]
    if (length(all_variants) == 0) return(GRanges())
  }

  # Annotate variants and filter down to a candidate set.
  benign_if_moot_statuses <- c(
    "Frameshifted",          # Background indel caused a frameshift.
    "Splice site disrupted", # Background variant broke a splice site.
    "DOWNSTREAM_OF_PTC"      # SNV falls after a newly introduced stop codon.
  )
  all_variants$CDS <- annotate_variants_with_cds(
    all_variants, variants_genomic_on_ts, txdb, genome,
    benign_if_moot_statuses)

  # A variant is considered "benign" if for ALL transcripts it affects,
  # the effect is either:
  #  1. Synonymous (no amino acid change).
  #  2. Moot (the transcript was already broken).
  #  3. Non-coding (the variant doesn't hit any CDS).
  is_valid_candidate <- sapply(all_variants$CDS, function(df) {
    # Case 1: Variant is non-coding (doesn't overlap any CDS). Keep it.
    if (nrow(df) == 0) {
      return(TRUE)
    }
    # Case 2: For coding variants, ALL its effects across all transcripts
    # must be benign.
    all(
      # Condition A: It's a synonymous mutation.
      # This will also catch NA == NA as OK
      (df$aa_ref == df$aa_alt) |
        # Condition B: Its effect is moot due to a compromised transcript.
        (df$status %in% benign_if_moot_statuses)
    )
  })
  all_variants <- all_variants[is_valid_candidate]
  if (length(all_variants) == 0) { return(GRanges()) }

  # is_known_variant column
  # 1. Best is synonymous codon that is compliant with genetic variant - true
  # 2. No genetic variant - NA
  # 3. Genetic variant non-compliant - false
  if (!is.null(snps)) {
    all_variants$dbSNP <- annotate_variants_with_snps(all_variants, snps)
  }
  all_variants$noncoding <- annotate_variants_with_noncoding(all_variants, annotation)
  if (!is.null(cadd)) {
    all_variants$CADD <- annotate_variants_with_cadd(all_variants, cadd)
  }

  if (alphagenome_key != "") {
    species <- if (organism(genome) == "Homo sapiens") {
      "human"
    } else if (organism(genome) == "Mus musculus") {
      "mouse"
    } else NA
    if (!is.na(species)) {
      ag_dt <- annotate_mutations_with_alphagenome(
        all_variants, alphagenome_key, species, python_exec, alphagenome_context)
      mcols(all_variants) <- cbind(mcols(all_variants), ag_dt)
    }
  }

  # Join annotated variants back to the guide map.
  # This step links each valid, annotated SNP with the specific guide(s)
  # it can disrupt.
  # We use `type = "equal"` to be precise about matching single-base locations.
  hits <- findOverlaps(all_variants, GRanges(candidate_snp_map$coords), type = "equal")
  if (length(hits) == 0) {
    # This could happen if, for example, all potential SNPs were filtered out
    return(GRanges())
  }

  # Construct the final data object.
  final_snps <- all_variants[S4Vectors::queryHits(hits)]
  final_snps$guide_name <- candidate_snp_map$guide_name[
    S4Vectors::subjectHits(hits)]
  final_snps$position_in_guide <- candidate_snp_map$position_in_guide[
    S4Vectors::subjectHits(hits)]
  names(final_snps) <- paste0(seqnames(final_snps), ":", start(final_snps), ":",
                              final_snps$REF, ">", final_snps$ALT,
                              " (for ", final_snps$guide_name, ")")

  all_statuses <- unlist(sapply(final_snps$CDS, `[[`, "status"))
  if (any(benign_if_moot_statuses %in% all_statuses)) {
    message(
      "NOTE: Your template includes indels that affect coding sequences. ",
      "Candidate SNPs have been generated in these regions and their predicted effects ",
      "are noted in the 'CDS' annotation column with statuses like 'frameshifted', ",
      "'downstream_of_ptc'. Please review these ",
      "candidates carefully, as standard synonymous criteria may not apply and ",
      "splicing may be affected."
    )
  }
  return(final_snps)
}
