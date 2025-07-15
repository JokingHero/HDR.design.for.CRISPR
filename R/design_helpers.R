#' Setup Design Environment
#' @description Loads and validates all initial transcript and genomic data.
#' @keywords internal
setup_design_environment <- function(
    ensemble_transcript_id, mutation_loci, mutation_original, annotation, genome) {
  # grab the transcript location on the genome
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(annotation))
  tx <- transcripts(txdb)
  tx <- tx[tx$tx_name == ensemble_transcript_id]
  cds <- get_cds(txdb, ensemble_transcript_id)
  mut_genomic <- get_genomic_mutation(cds, mutation_loci)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  stopifnot(paste0(getSeq(genome, cds[[1]]), collapse = "") == as.character(cds_seq))
  aa_cds_seq <- Biostrings::translate(cds_seq)

  if (isEmpty(mut_genomic)) {
    stop("Your `mutation_loci` does not seem to be on the CDS.")
  }
  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    stop("Your `mutation_original` is not the same as the one on the transcript sequence! Check transcirpt id.")
  }

  # genomic site is as big as whole gene transcript
  # we always show the tx sequence as 5'-3' (snapgene)
  # annotations also have to be flipped and be "+"
  genomic_seq <- Biostrings::getSeq(genome, tx) # this respects strand
  # for minus strand its reverseComplement of the forward - from gene perspective
  names(genomic_seq) <- ensemble_transcript_id

  # figure out mask - part of the template that actually is not part of the cds!!!
  # genome -------------------------------------
  # tx          ----------   ---- ---------
  # cds               ----   ---- ----
  # variant                    X
  cds_on_tx <- ranges(GenomicFeatures::pmapToTranscripts(
    cds[[1]], tx, ignore.strand = FALSE))
  stopifnot(cds_seq == genomic_seq[[1]][cds_on_tx])

  list(
    txdb = txdb, tx = tx, cds = cds, mut_genomic = mut_genomic,
    cds_seq = cds_seq, aa_cds_seq = aa_cds_seq, genomic_seq = genomic_seq,
    cds_on_tx = cds_on_tx)
}

#' Prepare Candidate Synonymous Mutations
#' @description Generates, filters, and annotates all possible synonymous SNPs.
#' @keywords internal
prepare_candidate_snps <- function(
    env_data, positions_to_mutate, mutation_loci, mutation_name,
    intron_bp, exon_bp, clinvar, snps, cadd, annotation, genome) {
  list2env(env_data, envir = environment())

  # ranges here are relative to the CDS
  mutations <- get_all_possilbe_mutations(positions_to_mutate, mutation_loci, cds_seq)
  # lets get genomic mutation positions
  mutations_genomic <- sapply(start(mutations), function(x) get_genomic_mutation(cds, x))
  mutations_genomic <- unlist(GRangesList(mutations_genomic))
  mutations <- GRanges(
    seqnames = mutation_name,
    ranges = ranges(GenomicFeatures::pmapToTranscripts( # tx positions
      mutations_genomic, tx, ignore.strand = FALSE)),
    strand = "+",
    original = elementMetadata(mutations)$original,
    replacement =  elementMetadata(mutations)$replacement,
    shift = elementMetadata(mutations)$shift,
    codon = elementMetadata(mutations)$codon,
    cds_pos = start(mutations)) # cds positions
  names(mutations) <- seq_along(mutations)
  if (!any(genomic_seq[[1]][ranges(mutations)] == DNAStringSet(mutations$original))) {
    stop("Something is not correct with generated synonymous-SNPs. Contact developers...")
  }

  # overlap mutations with anything on GFF
  is_not_over_ss <- over_splice_sites(
    mutations_genomic, txdb, intron_bp, exon_bp)
  mutations <- mutations[is_not_over_ss]
  mutations_genomic <- mutations_genomic[is_not_over_ss]
  if (!is.null(clinvar)) {
    vcf <- VariantAnnotation::readVcf(clinvar)
    vcf <- SummarizedExperiment::rowRanges(vcf)
    vcf <- vcf[lengths(vcf) == 1]
    seqlevelsStyle(vcf) <- seqlevelsStyle(mutations_genomic)
    is_not_over_clinvar <- !mutations_genomic %over% vcf
    mutations <- mutations[is_not_over_clinvar]
    mutations_genomic <- mutations_genomic[is_not_over_clinvar]
  }

  # compatible column
  # 1. Best is synonymous codon that is compliant with genetic variant - true
  # 2. No genetic variant - NA
  # 3. Genetic variant non-compliant - false
  if (!is.null(snps)) {
    mutations <- annotate_mutations_with_snps(
      mutations, mutations_genomic, snps)
  }
  mutations <- annotate_mutations_with_tx(genome, mutations, mutations_genomic, txdb)
  mutations <- annotate_mutations_with_noncoding(
    mutations, mutations_genomic, annotation)
  if (!is.null(cadd)) {
    mutations <- annotate_mutations_with_cadd(
      mutations, mutations_genomic, cadd)
  }

  # original mutation
  mut_tx <- ranges(GenomicFeatures::pmapToTranscripts( # tx positions
    mut_genomic, tx, ignore.strand = FALSE))
  names(mut_tx) <- mutation_name
  origin_mutation <- GRanges(
    seqnames = mutation_name,
    ranges = mut_tx,
    strand = "+",
    original = as.character(genomic_seq[[1]][mut_tx]),
    replacement = mutation_replacement, shift = 0,
    codon = ceiling(mutation_loci / 3))

  if (origin_mutation$original != mutation_original) {
    stop("Something is not ok with how our masks were calcualted here, report the error to the developers.")
  }

  return(list(mutations = mutations, mutations_genomic = mutations_genomic,
              mut_tx = mut_tx, origin_mutation = origin_mutation))
}

#' Primer 3
#' @description Design primers with the use of primer 3
#' @keywords internal
design_primers <- function(guide_distance,
         positions_to_mutate,
         genomic_seq,
         origin_mutation,
         mut_genomic,
         genome,
         primer3_path){
  # Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/full/path/to/primer3-x.x.x/src", sep = ":"))
  message("Designing primers...")

  max_dist <- max(guide_distance, abs(positions_to_mutate)) + 20
  options <- c(
    "PRIMER_SEQUENCE_ID",
    "SEQUENCE_TEMPLATE",
    "SEQUENCE_TARGET",
    "PRIMER_PICK_LEFT_PRIMER",
    "PRIMER_PICK_INTERNAL_OLIGO",
    "PRIMER_PICK_RIGHT_PRIMER",
    "PRIMER_OPT_SIZE",
    "PRIMER_MIN_SIZE",
    "PRIMER_MAX_SIZE",
    "PRIMER_MAX_NS_ACCEPTED",
    "PRIMER_PRODUCT_SIZE_RANGE",
    "P3_FILE_FLAG",
    "PRIMER_EXPLAIN_FLAG")
  values <- c(
    "template",
    as.character(genomic_seq[[1]]),
    paste0(start(origin_mutation) - max_dist - 1, ",", max_dist * 2),
    "1", "0", "1",
    getOption("PRIMER_OPT_SIZE", default = "22"),
    getOption("PRIMER_MIN_SIZE", default = "18"),
    getOption("PRIMER_MAX_SIZE", default = "25"),
    getOption("PRIMER_MAX_NS_ACCEPTED", default = "0"),
    paste0(getOption("PRODUCT_SIZE_MIN", default = "100"), "-",
           getOption("PRODUCT_SIZE_MAX", default = "250")),
    "0", "1"
  )
  opts <- c(paste0(options, "=", values), "=")

  out <- system2(command = primer3_path, input = opts,
                 stdout = TRUE)
  out <- strsplit(out, "=")
  out <- out[1:(length(out) - 1)]
  tags <- sapply(out, `[[`, 1)
  res <- sapply(out, `[[`, 2)

  pair_count <- as.numeric(res[tags == "PRIMER_PAIR_NUM_RETURNED"])
  primers <- list()
  for (i in seq_len(pair_count)) {
    sl <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_SEQUENCE")]
    sr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_SEQUENCE")]
    ll <- res[tags == paste0("PRIMER_LEFT_", i - 1)] #  "start,width"
    lr <- res[tags == paste0("PRIMER_RIGHT_", i - 1)]
    ll <- strsplit(ll, ",")[[1]]
    lr <- strsplit(lr, ",")[[1]]
    lw <- as.numeric(ll[2])
    rw <- as.numeric(lr[2])
    ll <- as.numeric(ll[1]) # coordinates already in BED format with -1
    lr <- as.numeric(lr[1])
    lr <- lr - rw + 1
    # left primer
    # genomic_seq[[1]][(ll + 1):(ll + lw)]
    # right primer
    # genomic_seq[[1]][(lr + 1):(lr + rw)]

    tmr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_TM")]
    tml <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_TM")]
    gcl <- res[tags == paste0("PRIMER_LEFT_", i - 1, "_GC_PERCENT")]
    gcr <- res[tags == paste0("PRIMER_RIGHT_", i - 1, "_GC_PERCENT")]
    ps <- res[tags == paste0("PRIMER_PAIR_", i - 1, "_PRODUCT_SIZE")]

    # offtargets detection with up to 2 mismatch
    ot_count <- 0
    for (ch in seqlevels(genome)) {
      # plus strand
      ot <- Biostrings::matchLRPatterns(
        sl, as.character(reverseComplement(DNAString(sr))),
        max.Lmismatch = 2, max.Rmismatch = 2,
        max.gaplength = 500, subject = genome[[ch]])
      # we need to filter too small products
      # we need to filter out repeated overlaps
      ot <- reduce(ranges(ot[nchar(ot) > (lw + rw)]))
      if (as.vector(seqnames(mut_genomic)) == ch) {
        ot <- setdiff(ot, IRanges(
          start(mut_genomic) - start(origin_mutation) + ll - lw,
          start(mut_genomic) - start(origin_mutation) + lr + rw + rw + 1))
      }
      if (length(ot) > 0) message(ch, "+")
      ot_count <- ot_count + length(ot)

      ot <- Biostrings::matchLRPatterns(
        sr, as.character(reverseComplement(DNAString(sl))),
        max.Lmismatch = 2, max.Rmismatch = 2,
        max.gaplength = 500, subject = genome[[ch]])
      ot <- reduce(ranges(ot[nchar(ot) > (lw + rw)]))
      if (as.vector(seqnames(mut_genomic)) == ch) {
        ot <- setdiff(ot, IRanges(
          start(mut_genomic) + start(origin_mutation) - lr - rw + 1 - rw,
          start(mut_genomic) + start(origin_mutation) - ll + lw))
      }
      if (length(ot) > 0) message(ch, "-")
      ot_count <- ot_count + length(ot)
    }

    primers[[i]] <- c(ps, ot_count, sl, ll, lw, tml, gcl, sr, lr, rw, tmr, gcr)
  }
  primers <- as.data.frame(do.call(rbind, primers))
  colnames(primers) <- c(
    "product_size", "offtarget_count",
    "seqeunce_left", "start_pos_left", "size_left",
    "melting_temp_left", "GC_left",
    "seqeunce_right", "start_pos_right", "size_right",
    "melting_temp_right", "GC_right")
  primers <- primers[order(primers$offtarget_count, primers$product_size), ]
  return(primers)
}

#' Export Design Results
#' @description Writes all final output files (FASTA, CSV, GFF3, GenBank).
#' @keywords internal
export_design_results <- function(
    output_dir, mutation_name, genomic_seq,
    origin_mutation, mutations, guides, repair_template,
    genome, cds_on_tx, aa_cds_seq,
    primers = NULL, probes_ = GRanges(), score_efficiency = FALSE,
    one_for_all = FALSE) {
  # prepare data to write
  # write bed/excel
  origin_mutation$codon <- NULL
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  all_combined <- c(origin_mutation, mutations, guides, repair_template, probes_)
  gxt <- data.frame("Guide name" = names(guides),
                    "Guide sequence" = as.vector(guides$original),
                    "Repair template name " = names(repair_template),
                    "Repair template sequence" = repair_template$replacement,
                    pam_disrupted = repair_template$pam_disrupted,
                    guide_disrupted = repair_template$guide_disrupted,
                    cadd = repair_template$cadd,
                    overlaps_something = repair_template$overlaps_something,
                    snp_quality = repair_template$snp_quality,
                    score_rank = if (score_efficiency) guides$rank_by_scores else NA)

  if (one_for_all) { # one template for all the guides
    gxt <- gxt[order(-gxt$score_rank, decreasing = T), ]
  } else { # each guide gets one template
    gxt <- gxt[order(!gxt$overlaps_something,
                     gxt$pam_disrupted, gxt$guide_disrupted, -gxt$cadd,
                     gxt$snp_quality,
                     -gxt$score_rank,
                     decreasing = T), ]
  }

  # Write FASTA file
  Biostrings::writeXStringSet(genomic_seq, file.path(output_dir, paste0(mutation_name, ".fa")))

  # Write BED/CSV file
  gr <- as.data.frame(all_combined)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  # Write guides x templates table
  write.table(
    gxt, file = file.path(output_dir, paste0(
      mutation_name, "_0based_guides_x_templates.csv")),
    quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(
    all_combined, file.path(output_dir, paste0(mutation_name, ".gff3")))
  write_genbank_custom(
    ensemble_transcript_id,
    genomic_seq[[1]],
    cds_on_tx,
    aa_cds_seq,
    organism(genome),
    file.path(output_dir, paste0(mutation_name, "_cds.gbk")))

  # Write primers if they exist
  if (!is.null(primers)) {
    write.table(
      primers,
      file = file.path(
        output_dir,
        paste0(mutation_name, "_primers.csv")),
      quote = F, sep = "\t", row.names = F, col.names = T)
  }
}
