# rm(list = ls(all.names = TRUE))
# gc(reset = T)
#
# library(Biostrings)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(IRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(crisprScore)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# library(VariantAnnotation)
#
# ensemble_transcript_id = "ENST00000399837.8"
# mutation_loci = 506
# mutation_original = "G"
# mutation_replacement = "A"
# mutation_name = "R169Q"
# output_dir = "~/"
# annotation = "/mnt/corsair/Projects/uib/CRIPSR/Oslo/gencode.v42.annotation.gff3.gz"
# genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# filter_to_guide = "TGCTGGAGGATTATCAGAAG"
# guide_distance = 10 + 17
# extension = 400
# positions_to_mutate = -19:19
# mutations_per_template = 0:3
# seed = 42
# score_efficiency = FALSE
# snps = SNPlocs.Hsapiens.dbSNP155.GRCh38
# clinvar = "/mnt/corsair/Projects/uib/CRIPSR/Oslo/clinvar.vcf.gz"
# intron_bp = 6
# exon_bp = 3
# probes = FALSE
# primer3 = "/home/ai/Soft/primer3/src/primer3_core" # or default ""
# source("./R/utils.R")

# ensemble_transcript_id = "ENST00000349496.11"
# mutation_loci = 110
# mutation_original = "C"
# mutation_replacement = "T"
# mutation_name = "C25118T"
#
# ensemble_transcript_id <- "ENST00000374202.7"
# mutation_loci <- 172
# mutation_original <- "C"
# mutation_replacement <- "T"
# mutation_name <- "pro58sel"
#
# # HAV2RC with the mutation c.245A>G (p.Tyr82Cys).
# ensemble_transcript_id <- "ENST00000307851.9"
# mutation_loci <- 245
# mutation_original <- "A"
# mutation_replacement <- "G"
# mutation_name <- "tyr82cys"
#
#
# ensemble_transcript_id = "ENST00000361099.8"
# mutation_loci = 820
# mutation_original = "C"
# mutation_replacement = "T"
# mutation_name = "STAT820"

#' @title Design experimental templates for each loci.
#' @description Designs all possible templates for each of the guides overlapping site of interest.
#' Templates are ranked based on: selecting synonymous SNPs that have the highest chance of disrupting
#' guides, minimizing disruption to other genomic functions by checking with submitted annotation file,
#' prioritizing SNPs that are already available in the population submitted through `snps
#' parameter. Guides are ranked first, by their inhibition potential, then by the quality of their snps,
#' and further by efficiency scoring.
#' @param ensemble_transcript_id This has to be ensemble transcript id e.g. ENST00000307851.9
#' @param mutation_loci Position on the gene that is mutated e.g. 245
#' @param mutation_original What was the original base on the genome? e.g. A
#' @param mutation_replacement What is the mutated base? e.g. G
#' @param mutation_name What is your name for this mutation e.g. tyr82cys
#' @param output_dir Where the files should be generated? Make sure you have writing permissions there.
#' @param annotation File path to the annotation file, a gff3.
#' @param genome BSgenome of your genome, compatible with annotation file, by default it is hg38.
#' @param filter_to_guide If "" then all gudies will be run, if you enter 20bp of your guide we will filter to only that guide.
#' @param guide_distance Window around which we should search for guides, relatively to the mutation loci. Defualt is 10 + 17bp.
#' @param extension How many bases upstream/downstream from the mutation loci should we include. Default is 400.
#' @param positions_to_mutate Which positions from the mutation are available for mutation. By default its -30:30
#' leaving 20bp on each side of the template for the in case of incomplete integration of the template. Also codon
#' occupied by the mutation is not included in the change.
#' @param mutations_per_template How many codons should be mutated per template. Insert a vector of
#' e.g. 0:5 which will then be used to generate all possible templates according to all allowed mutations per template combinations.
#' @param template_upstream How many bp upstream of mutation position you want to include in the repair template default: 59
#' @param template_downstream How many bp downstream of mutation position you want to include in the repair template default: 61
# '@param intron_bp How many bp of the intronic part of splicing should be excluded from synonymous SNPs (default 20, because of cryptic splice site potential)?
#' @param exon_bp How many bp of the exonic part of splicing should be excluded from synonymous SNPs (default 0)?
#' @param snps Can be either NULL or an object like SNPlocs.Hsapiens.dbSNP155.GRCh38 from library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
#' @param clinvar This is a clinVar database VCF file location. You can download it here https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz and https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
#' @param score_efficiency If you have installed relevant models from crisprScore package, you can set it to true to include scores from these models. Default is FALSE.
#' @param seed Ensures reproducibility of the random parts of the pipeline.
#' @param probes TRUE/FALSE whether to design probes or not?
#' @param primer3 If you install primer3 we can design primers, if you also input full path to primer3_core here.
#' @return writes files to the specified directory, might overwrite
#' @import Biostrings GenomicFeatures GenomicRanges SummarizedExperiment IRanges BSgenome BSgenome.Hsapiens.UCSC.hg38 VariantAnnotation GenomeInfoDb
#' @importFrom utils write.table combn
#' @export
#'
design_all_templates <- function(
    ensemble_transcript_id,
    mutation_loci,
    mutation_original,
    mutation_replacement,
    mutation_name,
    output_dir,
    annotation,
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    filter_to_guide = "",
    guide_distance = 10 + 17,
    extension = 400,
    positions_to_mutate = -19:19,
    mutations_per_template = 0:3,
    template_upstream = 59,
    template_downstream = 61,
    seed = 42,
    intron_bp = 6,
    exon_bp = 3,
    snps = NULL, # SNPlocs.Hsapiens.dbSNP155.GRCh38
    clinvar = NULL,
    score_efficiency = FALSE,
    probes = FALSE,
    primer3 = "") {
  set.seed(seed) # ensure reproducible randomness

  # grab the transcript location on the genome
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(annotation))
  cds <- get_cds(txdb, ensemble_transcript_id)
  mut_genomic <- get_genomic_mutation(cds, mutation_loci)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)[[1]]
  aa_cds_seq <- Biostrings::translate(cds_seq)

  if (isEmpty(mut_genomic)) {
    stop("Your `mutation_loci` does not seem to be on the CDS.")
  }
  if (as.character(Biostrings::getSeq(genome, mut_genomic)) != mutation_original) {
    stop("Your `mutation_original` is not the same as the one on the transcript sequence! Check transcirpt id.")
  }

  genomic_site <- resize(mut_genomic, width = extension * 2 + 1, fix = "center")
  genomic_seq <- Biostrings::getSeq(genome, genomic_site)
  names(genomic_seq) <- mutation_name

  # figure out mask - part of the template that actually is not part of the cds!!!
  mask_genome <- intersect(genomic_site, cds[[1]]) # this part on the genome is used for cds
  mask_cds <- ranges(GenomicFeatures::mapToTranscripts(mask_genome, cds)) # this part of the cds we are actually mutating
  # now - make it relative to the genome sequence we operate on
  tx_loci_with_introns <- resize(IRanges(mutation_loci + 1, width = 1),
                                 width = extension * 2, fix = "center")
  mask_seq <- shift(mask_cds, - (start(tx_loci_with_introns) - 1)) # parts of the sequence that are used for the cds

  # ranges here are relative to the extend
  mutations <- get_all_possilbe_mutations(
    mutation_name, positions_to_mutate, mutation_loci, cds_seq, aa_cds_seq, extension)

  # overlap mutations with anything on GFF
  # transform mutations to Genomic coordinates
  mutations_genomic <- sapply(mutations$cds_pos,
                              function(x) get_genomic_mutation(cds, x))
  mutations_genomic <- unlist(GRangesList(mutations_genomic))
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

  # original mutation
  origin_mutation <- GRanges(
    seqnames = mutation_name,
    ranges = IRanges(extension + 1, width = 1,
                     names = mutation_name),
    strand = "+",
    original = as.character(genomic_seq[[1]][extension + 1]),
    replacement = mutation_replacement, shift = 0,
    codon = ceiling(mutation_loci / 3))

  if (origin_mutation$original != mutation_original) {
    stop("Something is not ok with how our masks were calcualted here, report the error to the developers.")
  }

  # find and score guides in a window
  guides <- get_guides_and_scores(origin_mutation, mutation_name, guide_distance,
                                  genomic_seq, scores = score_efficiency)
  if (filter_to_guide != "") {
    found_guide <- guides$original == toupper(filter_to_guide)
    if (!any(found_guide)) {
      print(as.data.frame(guides))
      stop("Can't find this guide, see above for the guides we found.")
    }
    guides <- guides[found_guide]
  }
  pams <- resize(flank(guides, width = 3, start = FALSE), width = 2, fix = "end") # GG part
  probes_ <- GRanges()

  # now the repair template sequences
  repair_template <- GRanges()
  # prepare mutations ranges
  template_range <- promoters(ranges(origin_mutation),
                              upstream = template_upstream,
                              downstream = template_downstream)

  for (i in seq_along(guides)) { # each guide gets his own optimized mutations
    message("Working on guide: ", i)
    for (mpt in mutations_per_template) { # each number of mutations
      message("Mutations per template: ", mpt)
      if (mpt == 0) { # no mut
        mutated_seq <- genomic_seq[[1]]
        rt <- GRanges(seqnames = mutation_name,
                      ranges = template_range,
                      strand = "+",
                      replacement = as.character(extractAt(mutated_seq, template_range)),
                      pam_disrupted = 0,
                      guide_disrupted = 0,
                      overlaps_something = FALSE,
                      snp_quality = 0)
        temp_name <- paste0("Template_", as.character(names(guides[i])),
                            "_NoMut")
        names(rt) <- temp_name
        repair_template <- c(repair_template, rt)

        # design probes so that they overlap all of the HDR mutations
        if (probes) {
          candidates <- design_probes(mutation_name,
                                      start(origin_mutation) + min(positions_to_mutate),
                                      start(origin_mutation) + max(positions_to_mutate),
                                      mutated_seq,
                                      len_min = 19,
                                      len_max = 25,
                                      tmin = 50,
                                      tmax = 65,
                                      origin_mut_start = extension + 1)
          hdr_probes <- select_probes(origin_mutation, candidates, temp_name)
          names(hdr_probes) <- paste0("HDR probe ", seq_along(hdr_probes),
                                      " for ", temp_name)
          probes_ <- c(probes_, hdr_probes)
        }

      } else {
        muts <- get_all_combinations_of_mutations_for_guide(
          mutations, mpt, pams[i], guides[i])
        if (is.null(muts)) {
          message("Could not design synonymous SNPs for guide: ", guides[i], " ",
                  guides[i]$original," with ", mpt, " mutations per template")
          next
        }
        for (mut_i in muts) {
          mutated_seq <- replaceAt(genomic_seq[[1]],
                                   at = ranges(mut_i[[1]]),
                                   value = DNAStringSet(mut_i[[1]]$replacement))
          rt <- GRanges(seqnames = mutation_name,
                        ranges = template_range,
                        strand = "+",
                        replacement = as.character(extractAt(mutated_seq, template_range)),
                        pam_disrupted = mut_i[[2]],
                        guide_disrupted = mut_i[[3]],
                        overlaps_something = mut_i[[4]],
                        snp_quality = -1 * mut_i[[5]])
          temp_name <- paste0("Template_", as.character(names(guides[i])),
                              "_Mut_", paste0(names(mut_i[[1]]), collapse = "_"))
          names(rt) <- temp_name
          repair_template <- c(repair_template, rt)

          # design probes so that they overlap all of the HDR mutations
          if (probes) {
            candidates <- design_probes(mutation_name,
                                        start(origin_mutation) + min(positions_to_mutate),
                                        start(origin_mutation) + max(positions_to_mutate),
                                        mutated_seq,
                                        len_min = 19,
                                        len_max = 25,
                                        tmin = 50,
                                        tmax = 65,
                                        origin_mut_start = extension + 1)
            hdr_probes <- select_probes(c(mut_i[[1]], origin_mutation),
                                        candidates, temp_name)
            names(hdr_probes) <- paste0("HDR probe ", seq_along(hdr_probes),
                                        " for ", temp_name)
            probes_ <- c(probes_, hdr_probes)
          }
        }
      }
    }
  }

  if (probes) {
    # Reference probes - outside of the SNP area
    # Length/Tm= 20-25 bp, 60+/-1 oC
    ref_seq <- IRanges(start = start(origin_mutation), width = 1) + 120
    ref_seq <- setdiff(
      ref_seq, IRanges(start = start(origin_mutation), width = 1) +
        max(abs(positions_to_mutate)))
    # two chunks that can be used for Ref probes
    rp <- design_probes(
      mutation_name, start(ref_seq[1]), end(ref_seq[1]), genomic_seq[[1]],
      origin_mut_start = extension + 1)
    rp <- rp[which.max(rp$GC)]
    rp2 <- design_probes(
      mutation_name, start(ref_seq[2]), end(ref_seq[2]), genomic_seq[[1]],
      origin_mut_start = extension + 1)
    rp2 <- rp2[which.max(rp2$GC)]
    ref_probes <- c(rp, rp2)
    names(ref_probes) <- paste0("Ref probe ", seq_along(ref_probes))

    # NHEJ probe - ~20 bp and have Tm values 58-60oC
    # contain the original mutated sequence &
    # no SNPs and no correction of the mutation
    nhej_probes <- design_probes(
      mutation_name,
      start(origin_mutation)  + min(positions_to_mutate),
      start(origin_mutation)  + max(positions_to_mutate),
      genomic_seq[[1]],
      len_min = 19, len_max = 25, tmin = 50, tmax = 65,
      origin_mut_start = extension + 1)
    nhej_probes <- select_probes(origin_mutation, nhej_probes, temp_name)
    names(nhej_probes) <- paste0("NHEJ probe ", seq_along(nhej_probes))

    probes_ <- c(probes_, ref_probes, nhej_probes)
  }


  if (primer3 != "") {
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
      paste0(extension - max_dist - 1, ",", max_dist * 2),
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

    out <- system2(command = primer3, input = opts,
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
            start(mut_genomic) - extension + ll - lw,
            start(mut_genomic) - extension + lr + rw + rw + 1))
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
            start(mut_genomic) + extension - lr - rw + 1 - rw,
            start(mut_genomic) + extension - ll + lw))
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
    write.table(
      primers,
      file = file.path(
        output_dir,
        paste0(mutation_name, "_primers.csv")),
      quote = F, sep = "\t", row.names = F, col.names = T)
  }

  # first we write the sequence
  Biostrings::writeXStringSet(genomic_seq,
                              file.path(output_dir, paste0(mutation_name, ".fa")))
  # write bed/excel
  origin_mutation$codon <- NULL
  mutations$codon <- NULL
  # The only trick is remembering the BED uses 0-based coordinates. So add "-1" to the coords.
  all_combined <- c(origin_mutation, mutations, guides, repair_template)
  if (probes) {
    all_combined <- c(all_combined, probes_)
  }
  gr <- as.data.frame(all_combined)
  gr$start <- gr$start - 1
  write.table(gr, file = file.path(output_dir, paste0(mutation_name, "_0based.csv")),
              quote = F, sep = "\t", row.names = F, col.names = T)

  gxt <- data.frame("Guide name" = names(guides),
                    "Guide sequence" = guides$original,
                    "Repair template name " = names(repair_template),
                    "Repair template sequence" = repair_template$replacement,
                    pam_disrupted = repair_template$pam_disrupted,
                    guide_disrupted = repair_template$guide_disrupted,
                    overlaps_something = repair_template$overlaps_something,
                    snp_quality = repair_template$snp_quality,
                    score_rank = if (score_efficiency) guides$rank_by_scores else NA)

  gxt <- gxt[order(!gxt$overlaps_something,
                   gxt$pam_disrupted, gxt$guide_disrupted,
                   gxt$snp_quality,
                   -gxt$score_rank,
                   decreasing = T), ]
  write.table(
    gxt, file = file.path(output_dir, paste0(
      mutation_name, "_0based_guides_x_templates.csv")),
    quote = F, sep = "\t", row.names = F, col.names = T)

  # lets try to construct gff3 file for snapgene
  rtracklayer::export.gff3(
    all_combined, file.path(output_dir, paste0(mutation_name, ".gff3")))
  rtracklayer::export.gff3(
    GRanges(seqnames = mutation_name,
            ranges = mask_seq,
            strand = "+",
            type = "CDS"), file.path(output_dir, paste0(mutation_name, "_cds.gff3")))
}
