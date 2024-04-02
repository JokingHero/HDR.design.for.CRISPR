get_cds <- function(txdb, ensemble_transcript_id) {
  cds <- suppressWarnings(GenomicFeatures::cdsBy(txdb, by = "tx", use.names = T))
  cds <- cds[!duplicated(names(cds))]
  cds <- cds[ensemble_transcript_id]
  cds
}

get_genomic_mutation <- function(cds, mutation_loci) {
  mut_genomic <- GenomicFeatures::pmapFromTranscripts(
    IRanges(mutation_loci, width = 1), cds)
  mut_genomic <- mut_genomic[[1]][mut_genomic[[1]]$hit]
  mut_genomic
}

get_all_possilbe_mutations <- function(
    mutation_name,
    positions_to_mutate,
    mutation_loci,
    cds_seq,
    aa_cds_seq,
    extension) {
  # 3 extra mutations that are synonymous
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
  mutations <- GRanges()
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
    mut <- GRanges(seqnames = mutation_name,
                   ranges = IRanges(extension + i, width = 1),
                   strand = "+",
                   original = as.character(original_codon_seq[codon_position]),
                   replacement = "",
                   shift = i,
                   codon = codon_count,
                   cds_pos = this_pos_mutation_loci)
    mut <- rep(mut, length(replacement))
    mut$replacement <- replacement
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

get_combinations_of_mutations_for_guide <- function(mutations, mutations_per_template, pam, guide) {
  pam_disrupted <- which(ranges(mutations) %over% ranges(pam))
  guide_disrupted <- which(ranges(mutations) %over% ranges(guide))
  # not the same codons as in PAM
  if (length(pam_disrupted) != 0) {
    pam_disrupted <- sapply(
      unique(mutations$codon[pam_disrupted]),
      function(x) pam_disrupted[mutations$codon[pam_disrupted] == x][1])
  }
  if (length(guide_disrupted) != 0) {
    guide_disrupted <- guide_disrupted[
      !mutations$codon[guide_disrupted] %in% mutations$codon[pam_disrupted]]
    guide_disrupted <- sapply(
      unique(mutations$codon[guide_disrupted]),
      function(x) guide_disrupted[mutations$codon[guide_disrupted] == x][1])
  }

  if (length(pam_disrupted) >= mutations_per_template) {
    muts <- sample(pam_disrupted, mutations_per_template, replace = F)
  } else if ((length(pam_disrupted) + length(guide_disrupted)) >= mutations_per_template) {
    guide_disrupted <- guide_disrupted[
      order(distance(ranges(mutations[guide_disrupted]), ranges(pam)))]
    guide_disrupted <- guide_disrupted[1:(mutations_per_template - length(pam_disrupted))]
    muts <- c(pam_disrupted, guide_disrupted)
  } else {
    pam_and_guide <- c(pam_disrupted, guide_disrupted)
    npg <- seq_along(mutations)
    npg <- npg[!npg %in% pam_and_guide]
    npg <- sapply(
      unique(mutations$codon[npg]),
      function(x) npg[mutations$codon[npg] == x][1])
    npg <- npg[
      order(distance(ranges(mutations[npg]), ranges(pam)))]
    npg <- npg[1:(mutations_per_template - length(pam_and_guide))]
    muts <- c(pam_and_guide, npg)
  }
  list(muts, length(pam_disrupted), length(guide_disrupted))
}


annotate_mutations_with_snps <- function(mutations, mutations_genomic, snps) {
  seqlevelsStyle(mutations_genomic) <- seqlevelsStyle(snps)
  sbo <- GRanges(snpsByOverlaps(snps, mutations_genomic))
  hits <- findOverlaps(mutations_genomic, sbo)
  mutations$RefSNP_id <- ""
  mutations$alleles_as_ambig <- ""
  mutations$RefSNP_id[queryHits(hits)] <- sbo$RefSNP_id[subjectHits(hits)]
  mutations$alleles_as_ambig[queryHits(hits)] <- sbo$alleles_as_ambig[subjectHits(hits)]
  mutations$compatible <- mapply(function(base, iupac){
    stringr::str_detect(Biostrings:::IUPAC_CODE_MAP[iupac], base)
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

annotate_mutations_with_tx <- function(mutations, mutations_genomic, txdb) {
  mutations$syn_tx_count <- 0
  mutations$syn_tx <- ""
  mutations$overlap_tx_count <- 0

  all_tx_cds <- suppressWarnings(GenomicFeatures::cdsBy(
    txdb, by = "tx", use.names = T))
  for (i in seq_along(mutations)) {
    to_check <- all_tx_cds[all_tx_cds %over% mutations_genomic[i]]
    i_cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, to_check)
    i_aa_cds_seq <- Biostrings::translate(i_cds_seq)
    i_loc <- GenomicFeatures::mapToTranscripts(mutations_genomic[i], to_check)
    nonsyn <- c()
    for (j in seq_along(to_check)) {
      ji_cds_seq <- i_cds_seq[[j]]
      if (as.character(ji_cds_seq[start(i_loc[j])]) != mutations[i]$original) {
        stop("Nonysynomous SNPs reference mismatch.")
      }
      ji_cds_seq[start(i_loc[j])] <- DNAString(mutations[i]$replacement)
      if (i_aa_cds_seq[[j]] != Biostrings::translate(ji_cds_seq)) {
        nonsyn <- c(nonsyn, names(to_check)[j])
      }
    }
    mutations$syn_tx_count[i] <- length(nonsyn)
    mutations$syn_tx[i] <- paste0(nonsyn, collapse = ";")
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

melting_temp <- function(x) TmCalculator::Tm_GC(
  ntseq = as.character(x),
  variant = "Primer3Plus", Na = 50, outlist = FALSE)

gc_fract <- function(x) letterFrequency(x, letters = "CG", as.prob = TRUE)

design_probes <- function(mutation_name, st, sp, s, tmin = 59, tmax = 61, len_min = 20, len_max = 25,
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
