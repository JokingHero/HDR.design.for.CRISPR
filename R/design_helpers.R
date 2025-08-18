#' Primer 3
#' @description Design primers with the use of primer 3
#' We want at least one primer to be outside of the template_range
#' We use `template_range_extended` as the region where we can design
#' We check for off-targets
#' We want to make sure the target_on_template_extended is captured inside the PCRed sequence
#' @keywords internal
design_primers <- function(
    template_range_extended, template_range, variant_with_allowed, genome, primer3_path){
  message("Designing primers...")
  options <- c(
    "PRIMER_SEQUENCE_ID",
    "SEQUENCE_TEMPLATE",
    "SEQUENCE_TARGET",
    "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST",
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
  template_range_extended_seq <- getSeq(genome, template_range_extended)[[1]]
  variant_with_allowed_template_extended <- pmapToTranscripts(
    variant_with_allowed, template_range_extended)
  template_range_on_template_extended <- pmapToTranscripts(
    template_range, template_range_extended)
  values <- c(
    "template",
    as.character(template_range_extended_seq),
    paste0(start(variant_with_allowed_template_extended), ",",
           width(variant_with_allowed_template_extended)),
    paste0(1, ",", start(template_range_on_template_extended), ",, ; ,,",
           end(template_range_on_template_extended), ",",
           nchar(template_range_extended_seq) - end(template_range_on_template_extended) - 1),
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
    # template_range_extended_seq[(ll + 1):(ll + lw)]
    # right primer
    # template_range_extended_seq[(lr + 1):(lr + rw)]

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
      if (as.vector(seqnames(template_range_extended)) == ch) {
        ot <- setdiff(ot, ranges(template_range_extended))
      }
      ot_count <- ot_count + length(ot)

      ot <- Biostrings::matchLRPatterns(
        sr, as.character(reverseComplement(DNAString(sl))),
        max.Lmismatch = 2, max.Rmismatch = 2,
        max.gaplength = 500, subject = genome[[ch]])
      ot <- reduce(ranges(ot[nchar(ot) > (lw + rw)]))
      if (as.vector(seqnames(template_range_extended)) == ch) {
        ot <- setdiff(ot, ranges(template_range_extended))
      }
      if (length(ot) > 0) message(ch, "-")
      ot_count <- ot_count + length(ot)
    }

    primers[[i]] <- c(ps, ot_count, sl, ll, lw, tml, gcl, sr, lr, rw, tmr, gcr)
  }
  primers <- as.data.frame(do.call(rbind, primers))
  colnames(primers) <- c(
    "product_size", "offtarget_count",
    "LEFT_sequence", "LEFT_start_pos", "LEFT_size",
    "LEFT_melting_temp", "LEFT_GC",
    "RIGHT_sequence", "RIGHT_start_pos", "RIGHT_size",
    "RIGHT_melting_temp", "RIGHT_GC")
  primers <- primers[order(primers$offtarget_count, primers$product_size), ]
  primers_genomic_l <- pmapFromTranscripts(IRanges(start = as.numeric(primers$LEFT_start_pos) + 1,
                                                   width = as.numeric(primers$LEFT_size)),
                                           template_range_extended)
  primers_genomic_r <- pmapFromTranscripts(IRanges(start = as.numeric(primers$RIGHT_start_pos) + 1,
                                                   width = as.numeric(primers$RIGHT_size)),
                                           template_range_extended)
  primers$LEFT_genomic_start <- start(primers_genomic_l)
  primers$RIGHT_genomic_start <- start(primers_genomic_r)
  primers$chrom <- as.character(seqnames(primers_genomic_l))
  primers$rownames <- paste0("PAIR ", seq_along(primers_genomic_l))
  return(primers)
}

#' @title Helper write components
#' @description Helper to write CSV and GFF3 for a given GRanges component.
#' It calculates relative coordinates within the `chrom_relative` window.
#' @keywords internal
#'
write_component_files <- function(
    output_dir, design_name, component, component_name, chrom_relative) {
  if (length(component) == 0) {
    return(invisible(NULL))
  }

  # Calculate relative coordinates within the chrom_relative window
  component_relative <- pmapToTranscripts(component, chrom_relative)
  mcols(component)$start_relative <- start(component_relative)
  mcols(component)$end_relative <- end(component_relative)

  # Write 1-based coordinate CSV file
  df_component <- as.data.frame(component)
  df_component$row_names <- names(component)
  write.table(
    df_component,
    file = file.path(output_dir, paste0(design_name, "_", component_name, "_1based.csv")),
    quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  rtracklayer::export.gff3(
    component,
    file.path(output_dir, paste0(design_name, "_", component_name, ".gff3")))
}

#' Export Design Results
#' @description Writes all final output files (FASTA, CSV, GFF3, GenBank).
#' @keywords internal
export_design_results <- function(
    output_dir, design_name,
    variant_genomic,
    var_data, guides, repair_template,
    genome, txdb,
    primers = GRanges(),
    probes_ = GRanges(),
    score_efficiency = FALSE,
    one_for_all = FALSE) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # create large region of 1k BP arond the variant and map all coordinates there
  chrom_relative <- promoters(
    variant_genomic, 1000, 1000 + width(variant_genomic), use.names = FALSE)
  chrom_relative_seq <- getSeq(genome, chrom_relative)[[1]]
  strand(chrom_relative) <- "*"
  # Write FASTA file
  chrom <- DNAStringSet(chrom_relative_seq)
  names(chrom) <- paste0(chrom_relative)
  Biostrings::writeXStringSet(
    DNAStringSet(chrom_relative_seq), file.path(output_dir, paste0(design_name, ".fa")))

  if (length(guides) > 0) {
    mcols(guides) <- mcols(guides)[
      , !apply(mcols(guides), 2, function(x) all(is.na(x))), drop = FALSE]
  }

  # Write files for each component using the helper function
  write_component_files(output_dir, design_name, variant_genomic, "query", chrom_relative)
  write_component_files(output_dir, design_name, guides, "guides", chrom_relative)
  repair_template <- repair_template[
    order(!repair_template$overlaps_noncoding,
          repair_template$pam_disrupted,
          repair_template$guide_disrupted,
          -repair_template$cadd,
          repair_template$snp_quality,
          decreasing = T),
  ]
  write_component_files(output_dir, design_name, repair_template, "templates", chrom_relative)
  write_component_files(output_dir, design_name, probes_, "probes", chrom_relative)

  if (!is.null(primers)) {
    # Calculate relative coordinates within the chrom_relative window
    primers_left <- GRanges(seqnames = primers$chrom,
                            ranges = IRanges(start = as.numeric(primers$LEFT_genomic_start),
                                             width = as.numeric(primers$LEFT_size)),
                            strand = "+")
    left_relative <- pmapToTranscripts(primers_left, chrom_relative)
    primers$LEFT_start_pos <- start(left_relative)
    primers_right <- GRanges(seqnames = primers$chrom,
                            ranges = IRanges(start = as.numeric(primers$RIGHT_genomic_start),
                                             width = as.numeric(primers$RIGHT_size)),
                            strand = "+")
    right_relative <- pmapToTranscripts(primers_right, chrom_relative)
    primers$RIGHT_start_pos <- start(right_relative)
    write.table(
      as.data.frame(primers),
      file = file.path(output_dir, paste0(design_name, "_primers_1based.csv")),
      quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
  }

  # We move the complex annotations to separte tables
  var_data_cds <- as.data.frame(var_data$CDS)
  if (nrow(var_data_cds) > 0) var_data_cds$group_name <- names(var_data)[var_data_cds$group]
  var_data_dbSNP <- as.data.frame(var_data$dbSNP)
  if (nrow(var_data_dbSNP) > 0) var_data_dbSNP$group_name <- names(var_data)[var_data_dbSNP$group]
  var_data_noncoding <- as.data.frame(var_data$noncoding)
  if (nrow(var_data_noncoding) > 0) var_data_noncoding$group_name <-
    names(var_data)[var_data_noncoding$group]
  var_data$CDS <- NULL
  var_data$dbSNP <- NULL
  var_data$noncoding <- NULL
  write_component_files(output_dir, design_name, var_data, "variants", chrom_relative)

  if (!isEmpty(var_data_cds)) {
    write.table(var_data_cds, file = file.path(
      output_dir, paste0(design_name, "_cds_annotations_1based.csv")),
      quote = F, sep = ",", row.names = F, col.names = T)

    cds <- cdsBy(txdb, by = "tx", use.names = T)
    cds <- cds[names(cds) %in% var_data_cds$tx_id]
    cds <- restrict(cds, start = start(chrom_relative), end = end(chrom_relative), use.names = T)
    cds <- unlist(cds, use.names = T)
    cds_relative <- pmapToTranscripts(cds, chrom_relative)
    cds$start_relative <- start(cds_relative)
    cds$end_relative <- end(cds_relative)
    cds$row_names <- names(cds)
    names(cds) <- NULL
    cds_dt <- as.data.frame(cds)
    write.table(cds_dt, file = file.path(
      output_dir, paste0(design_name, "_cds_1based.csv")),
      quote = F, sep = ",", row.names = F, col.names = T)

    write_genbank_custom(
      paste0(chrom_relative),
      chrom_relative_seq,
      cds,
      NULL, # TODO fix this amino acids annotation
      organism(genome),
      file.path(output_dir, paste0(design_name, "_cds.gbk")))
  }

  if (!isEmpty(var_data_dbSNP)) {
    write.table(var_data_dbSNP, file = file.path(
      output_dir, paste0(design_name, "_dbSNP_annotations_1based.csv")),
      quote = F, sep = ",", row.names = F, col.names = T)
  }

  if (!isEmpty(var_data_noncoding)) {
    write.table(var_data_noncoding, file = file.path(
      output_dir, paste0(design_name, "_noncoding_annotations_1based.csv")),
      quote = F, sep = ",", row.names = F, col.names = T)
  }
}
