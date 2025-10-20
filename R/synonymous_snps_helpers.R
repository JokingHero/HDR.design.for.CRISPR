#' Define Genomic Search Positions Relative to a Variant
#'
#' @description
#' This helper function calculates the absolute genomic coordinates for a set of
#' positions defined relatively (upstream and downstream) to a given variant.
#'
#' The function operates on the reference genome's coordinate system. For
#' deletions, "upstream" is relative to the start of the deleted block and
#' "downstream" is relative to the end. For insertions, both are relative to
#' the single anchor base on the reference. The strand of the variant is ignored,
#' as positions are calculated on the genomic sequence.
#'
#' @param variant_genomic A GRanges object of length 1 representing the
#'   central variant. The start and end coordinates are used as anchors.
#' @param allowed_positions_upstream A numeric vector of positive integers
#'   indicating the distances to search upstream from the variant's start
#'   position (e.g., `c(1, 5)` means 1bp and 5bp upstream).
#' @param allowed_positions_downstream A numeric vector of positive integers
#'   indicating the distances to search downstream from the variant's end
#'   position (e.g., `c(2, 10)` means 2bp and 10bp downstream).
#'
#' @return A sorted GRanges object containing all the specified genomic
#'   positions, each with a width of 1. Returns an empty GRanges object if
#'   no positions are specified.
#' @keywords internal
#'
get_search_positions <- function(variant_genomic,
                                 allowed_positions_upstream,
                                 allowed_positions_downstream) {
  if (!methods::is(variant_genomic, "GRanges") || length(variant_genomic) != 1) {
    stop("`variant_genomic` must be a GRanges object of length 1.", call. = FALSE)
  }
  if (any(allowed_positions_upstream <= 0)) {
    stop("`allowed_positions_upstream` must be a vector of positive numbers.", call. = FALSE)
  }
  if (any(allowed_positions_downstream <= 0)) {
    stop("`allowed_positions_downstream` must be a vector of positive numbers.", call. = FALSE)
  }

  # Upstream positions are relative to the start of the variant on the reference
  upstream_coords <- GenomicRanges::start(variant_genomic) - allowed_positions_upstream
  # Downstream positions are relative to the end of the variant on the reference
  downstream_coords <- GenomicRanges::end(variant_genomic) + allowed_positions_downstream
  search_coords <- c(upstream_coords, downstream_coords)

  if (length(search_coords) == 0) {
    return(GenomicRanges::GRanges())
  }
  potential_snps <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(variant_genomic),
    ranges = IRanges::IRanges(start = search_coords, width = 1),
    strand = "*")
  return(sort(potential_snps))
}

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

#' Annotate coding sequence changes for a set of variants.
#'
#' @description
#' This function takes a GRanges object of single nucleotide variants (SNVs),
#' all on the '+' strand, and annotates their effects on the coding sequences (CDS)
#' of overlapping transcripts.
#' @param all_variants A GRanges object containing the variants. Must have "REF"
#'   and "ALT" metadata columns and be named.
#' @param genome A BSgenome object for the reference genome.
#' @param txdb A TxDb object containing transcript models.
#' @return A DataFrameList where each element corresponds to a variant. Each
#'   DataFrame contains detailed annotations for every transcript the variant
#'   overlaps, with columns: tx_id, tx_strand, codon_ref, codon_alt,
#'   aa_ref, aa_alt, codon_num, and frame.
#' @keywords internal
#'
annotate_variants_with_cds <- function(all_variants, txdb, genome) {
  cds_by_tx <- suppressWarnings(GenomicFeatures::cdsBy(
    txdb, by = "tx", use.names = TRUE
  ))
  cds_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, cds_by_tx)
  result <- lapply(seq_along(all_variants), function(i) {
    variant <- all_variants[i]
    strand(variant) <- "*"
    hits <- suppressWarnings(
      GenomicFeatures::mapToTranscripts(variant, cds_by_tx, ignore.strand = FALSE))
    if (length(hits) == 0) {
      return(S4Vectors::DataFrame())
    }

    # Get transcript IDs and their actual strands
    tx_ids <- as.character(seqnames(hits))
    tx_strands <- as.character(strand(hits))

    is_minus_strand <- tx_strands == "-"
    effective_alt <- rep(variant$ALT, length(hits))
    if (any(is_minus_strand)) {
      rc_alt <- as.character(reverseComplement(DNAString(variant$ALT)))
      effective_alt[is_minus_strand] <- rc_alt
    }
    # Get position, frame, and codon number within the transcript's CDS
    pos_in_cds <- start(hits)
    frame <- (pos_in_cds - 1) %% 3 + 1
    codon_num <- floor((pos_in_cds - 1) / 3) + 1

    # Extract the reference codon for every hit
    ref_codons_dna <- mapply(function(tid, cnum) {
      transcript_seq <- cds_seqs[[tid]]
      codon_start_pos <- (cnum - 1) * 3 + 1
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
    }, tx_ids, codon_num, SIMPLIFY = FALSE)
    ref_codons_dna <- DNAStringSet(ref_codons_dna)
    alt_codons_dna <- replaceAt(
      ref_codons_dna,
      IRangesList(lapply(frame, IRanges, width = 1)),
      CharacterList(split(effective_alt, seq_along(effective_alt))))

    S4Vectors::DataFrame(
      tx_id = tx_ids,
      ranges = start(hits),
      tx_strand = tx_strands,
      codon_ref = as.character(ref_codons_dna),
      codon_alt = as.character(alt_codons_dna),
      aa_ref = translate_safe(ref_codons_dna),
      aa_alt = translate_safe(alt_codons_dna),
      codon_num = codon_num,
      frame = frame)
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
                                                python_exec = "python3") {
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
  mean_raw_score <- raw_score <- variant_id <- output_type <- NULL
  ag_dt <- readr::read_csv(temp_output_csv, show_col_types = FALSE)
  # we will now average all raw scores for many tissues and cell lines
  ag_dt <- dplyr::group_by(ag_dt, variant_id, output_type)
  ag_dt <- suppressMessages(dplyr::summarise(
    ag_dt, mean_raw_score = mean(raw_score, na.omit = T)))
  ag_dt$output_type <- paste0("alphagenome_", ag_dt$output_type)
  ag_dt <- tidyr::pivot_wider(
    ag_dt,
    id_cols = variant_id,
    names_from = output_type,
    values_from = mean_raw_score)

  original_ids_order <- paste0(
    as.character(GenomicRanges::seqnames(all_variants)), ":",
    GenomicRanges::start(all_variants), ":",
    all_variants$REF, ">",
    all_variants$ALT)

  ag_dt <- ag_dt[match(ag_dt$variant_id, original_ids_order), ]
  ag_dt$variant_id <- NULL
  return(ag_dt)
}

#' Prepare Candidate Synonymous Mutations
#' @description Generates, filters, and annotates all possible synonymous SNPs.
#' @keywords internal
prepare_candidate_snps <- function(
    annotation, txdb, genome,
    variant_genomic, allowed_positions_upstream, allowed_positions_downstream,
    intron_bp, exon_bp, clinvar, snps, cadd,
    alphagenome_key, python_exec) {

  search_positions <- get_search_positions(
    variant_genomic = variant_genomic,
    allowed_positions_upstream = allowed_positions_upstream,
    allowed_positions_downstream = allowed_positions_downstream)
  if (length(search_positions) == 0) stop("Unable to design for selected positions.", call. = F)

  all_variants <- get_all_possible_mutations(search_positions, genome)
  all_variants$CDS <- annotate_variants_with_cds(all_variants, txdb, genome)

  # Only benign SNPs are left
  # empty DataFrame for no overlapping CDS also is a pass
  all_variants <- all_variants[sapply(all_variants$CDS, function(x) all(x$aa_ref == x$aa_alt))]
  # filter out too close to the splice sites
  is_outside <- is_outside_splice_sites(all_variants, txdb, intron_bp, exon_bp)
  all_variants <- all_variants[is_outside]

  if (!is.null(clinvar)) {
    vcf <- VariantAnnotation::readVcf(clinvar)
    vcf <- SummarizedExperiment::rowRanges(vcf)
    vcf <- vcf[lengths(vcf) == 1]
    seqlevelsStyle(vcf) <- seqlevelsStyle(all_variants)
    is_not_over_clinvar <- !(all_variants %over% vcf)
    all_variants <- all_variants[is_not_over_clinvar]
  }

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
        all_variants, alphagenome_key, species, python_exec)
      mcols(all_variants) <- cbind(mcols(all_variants), ag_dt)
    }
  }

  names(all_variants) <- paste0(seqnames(all_variants), ":", start(all_variants), ":",
                                all_variants$REF, ">", all_variants$ALT)
  return(all_variants)
}
