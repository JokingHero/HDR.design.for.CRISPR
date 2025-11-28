#' @keywords internal
#'
parse_primer3_system2_output <- function(output_lines) {

  # 1. Identify where records end (lines that are just "=")
  # We use trimws to ensure we catch it even if there are spaces like " ="
  sep_indices <- which(trimws(output_lines) == "=")

  if (length(sep_indices) == 0) {
    warning("No '=' separators found in output. Primer3 might have failed or output format is wrong.")
    return(data.frame())
  }

  # 2. Define Start and End indices for each record
  # First record starts at 1. Subsequent records start after the previous '='.
  start_indices <- c(1, head(sep_indices, -1) + 1)
  end_indices   <- sep_indices - 1 # Exclude the '=' line itself

  all_pairs_list <- list()

  # 3. Iterate through records
  for (i in seq_along(start_indices)) {
    # Extract the subset of lines for this single record
    # Ensure we don't index weirdly if start > end (empty record)
    if (start_indices[i] > end_indices[i]) next

    rec_lines <- output_lines[start_indices[i]:end_indices[i]]

    # Filter only lines containing "=" to avoid parsing garbage/warnings
    rec_lines <- rec_lines[grepl("=", rec_lines)]
    if (length(rec_lines) == 0) next

    # 4. Fast Vectorized Tag Parsing
    # Find position of first '='
    eq_pos <- regexpr("=", rec_lines, fixed = TRUE)

    keys <- substr(rec_lines, 1, eq_pos - 1)
    vals <- substr(rec_lines, eq_pos + 1, nchar(rec_lines))
    tags <- stats::setNames(vals, keys)

    # 5. Extract Data
    n_ret <- as.numeric(tags["PRIMER_PAIR_NUM_RETURNED"])

    if (!is.na(n_ret) && n_ret > 0) {
      seq_id    <- tags["SEQUENCE_ID"]
      int_oligo <- tags["SEQUENCE_INTERNAL_OLIGO"]

      # Iterate over pairs (0 to N-1)
      for (k in 0:(n_ret - 1)) {

        # Helper to safely get tag
        get_val <- function(key) tags[paste0(key, "_", k)]
        get_sub_val <- function(key, suffix) tags[paste0(key, "_", k, suffix)]

        # Parse Coordinates "start,length"
        left_coords  <- strsplit(get_val("PRIMER_LEFT"), ",")[[1]]
        right_coords <- strsplit(get_val("PRIMER_RIGHT"), ",")[[1]]

        # Build Dataframe Row
        pair_df <- data.frame(
          SEQUENCE_ID             = if(is.na(seq_id)) NA else seq_id,
          SEQUENCE_INTERNAL_OLIGO = if(is.na(int_oligo)) NA else int_oligo,
          PRIMER_PAIR_NUMBER      = k,

          PRIMER_FORWARD_SEQUENCE    = get_sub_val("PRIMER_LEFT", "_SEQUENCE"),
          PRIMER_REVERSE_SEQUENCE   = get_sub_val("PRIMER_RIGHT", "_SEQUENCE"),

          # Convert to Numeric
          PRIMER_FORWARD_START       = as.numeric(left_coords[1]),
          PRIMER_REVERSE_START      = as.numeric(right_coords[1]),
          PRIMER_FORWARD_LENGTH      = as.numeric(left_coords[2]),
          PRIMER_REVERSE_LENGTH     = as.numeric(right_coords[2]),

          PRIMER_FORWARD_TM          = as.numeric(get_sub_val("PRIMER_LEFT", "_TM")),
          PRIMER_REVERSE_TM         = as.numeric(get_sub_val("PRIMER_RIGHT", "_TM")),

          PRODUCT_SIZE            = as.numeric(get_sub_val("PRIMER_PAIR", "_PRODUCT_SIZE")),

          PRIMER_FORWARD_GC          = as.numeric(get_sub_val("PRIMER_LEFT", "_GC_PERCENT")),
          PRIMER_REVERSE_GC         = as.numeric(get_sub_val("PRIMER_RIGHT", "_GC_PERCENT")),

          stringsAsFactors = FALSE
        )

        all_pairs_list[[length(all_pairs_list) + 1]] <- pair_df
      }
    }
  }

  # 6. Final Combine
  if (length(all_pairs_list) == 0) {
    return(data.frame())
  }

  return(do.call(rbind, all_pairs_list))
}

#' Primer 3
#' @description Design primers with the use of primer 3 -
#' We want at least one primer to be outside of the template_range
#' We use `template_range_extended` as the region where we can design
#' We check for off-targets
#' We want to make sure the target_on_template_extended is captured inside the PCRed sequence
#' @keywords internal
#'
design_primers_rna <- function(
    target,
    position,
    guides_seq,
    primers_per_guide,
    primer3_path){
  message("Designing primers...")
  global_str <- list(
    "PRIMER_TASK" = "generic",
    "PRIMER_FIRST_BASE_INDEX" = 1, # Very important for R - primer3 compatibility
    "PRIMER_PICK_LEFT_PRIMER" = 1,
    "PRIMER_PICK_RIGHT_PRIMER" = 1,
    "PRIMER_PICK_INTERNAL_OLIGO" = 1, # We are treating the Guide as the Internal Oligo
    # Guide Constraints (Fixed at 28)
    "PRIMER_INTERNAL_MIN_SIZE" = 28,
    "PRIMER_INTERNAL_OPT_SIZE" = 28,
    "PRIMER_INTERNAL_MAX_SIZE" = 28,
    # Force Primer3 to accept the guide even if Tm is weird
    "PRIMER_PICK_ANYWAY" = 1,
    # Relax internal oligo constraints so P3 doesn't fail just because of the guide
    "PRIMER_INTERNAL_MIN_TM" = 0,
    "PRIMER_INTERNAL_MAX_TM" = 100,
    "PRIMER_INTERNAL_MIN_GC" = 0,
    "PRIMER_NUM_RETURN" = primers_per_guide,
    "PRIMER_PRODUCT_SIZE_RANGE" = paste0(getOption("PRODUCT_SIZE_MIN", default = "75"), "-",
                                         getOption("PRODUCT_SIZE_MAX", default = "150")),
    "PRIMER_MIN_SIZE" = getOption("PRIMER_MIN_SIZE", default = "18"),
    "PRIMER_OPT_SIZE" = getOption("PRIMER_OPT_SIZE", default = "20"),
    "PRIMER_MAX_SIZE" = getOption("PRIMER_MAX_SIZE", default = "27"),
    "PRIMER_MIN_TM" = 55,
    "PRIMER_OPT_TM" = 60,
    "PRIMER_MAX_TM" = 65,
    "P3_FILE_FLAG" = 0)

  input_records <- character(length(guides_seq))
  for (i in seq_along(guides_seq)) {
    # SEQUENCE_INTERNAL_OLIGO forces P3 to use THIS sequence as the probe
    input_records[i] <- paste0(
      "SEQUENCE_ID=", i, "\n",
      "SEQUENCE_TEMPLATE=", target, "\n",
      "SEQUENCE_INTERNAL_OLIGO=", as.character(guides_seq[[i]]), "\n",
      if (i == length(guides_seq)) "=" else "=\n")
  }

  # BoulderIO allows Globals to persist.
  final_input_str <- c(
    paste0(paste0(names(global_str), "=", global_str, collapse = "\n"),
           "\n",
           paste(input_records, collapse = "")))

  out <- tryCatch({
    system2(command = primer3_path,
            input = final_input_str, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    stop("Failed to execute primer3. Error: ", e$message)
  })

  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("Primer3 execution failed with exit code ", status, ". Output:\n", paste(out, collapse="\n"))
  }

  if (is.null(out) || length(out) < 2) {
    stop("Primer3 did not return valid output.")
  }

  primers <- parse_primer3_system2_output(out)
  rownames(primers) <- paste0("Design ", 1:nrow(primers))
  # TODO off-target search for primers?! can we even do that?!
  return(primers)
}

#' Predict Cas13a Guide-Target Interaction using RNAup
#'
#' @description Calculates the thermodynamic properties of the interaction
#' between Cas13a guides and PCR amplicons. It considers both the hybridization
#' energy and the accessibility of the target site within the amplicon.
#'
#' @param spacer_with_linker_seq A RNAStringSet or character vector of spacer sequences with LINKER!
#' @param amp_seq_rna A RNAStringSet of amplicon sequences
#' @param rnaup_path Path to the RNAup binary.
#' @param window_size Integer. Max interaction length (roughly spacer length). Default 30.
#'
#' @return A data.frame containing:
#' \item{dG_binding}{Total interaction energy (lower is better)}
#' \item{dG_duplex}{Energy of hybridization}
#' \item{dG_opening}{Energy required to open the target structure}
#' @keywords internal
#'
predict_spacer_amplicon_interaction <- function(
    spacer_with_linker_seq, amp_seq_rna,
    rnaup_path = "RNAup",
    window_size = 30) {
  if (length(spacer_with_linker_seq) != length(amp_seq_rna)) {
    stop("Spacers or Amplicons must be 1 each.")
  }

  # Prepare input for RNAup
  # RNAup with -b -w
  # Format: >ID \n SeqLong&SeqShort
  ids <- paste0(seq_along(spacer_with_linker_seq))
  input_lines <- paste0(">", ids, "\n", amp_seq_rna, "&", spacer_with_linker_seq)
  # -b: include probability of being unpaired for both RNA
  # --noLP: No lonely pairs (often gives more realistic biophysical results)
  args <- c("-b", "-w", as.character(window_size), "--no_output_file", "--noLP")

  results_list <- list()
  for (i in 1:length(amp_seq_rna)) {
    out <- tryCatch({
      system2(command = rnaup_path, args = args,
              input = input_lines[i], stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      out <- ""
    })

    # Look for the line containing " = "
    data_line_idx <- grep(" = ", out, fixed = TRUE)

    if (length(data_line_idx) > 0) {
      raw_line <- out[data_line_idx]

      # --- PARSING STRATEGY (Split by Delimiters) ---

      # Structure of line:
      # "((..))  Start,End  :  Start,End  (dG_Bind = dG_Duplex + dG_Open1 + dG_Open2)"

      # Step A: Split into Coords Part (Left) and Energy Part (Right) using " ("
      parts_main <- strsplit(raw_line, " (", fixed = TRUE)[[1]]
      coords_part <- parts_main[1]  # "((..))  Start,End  :  Start,End"
      energy_part <- parts_main[2]  # "-30.10 = -39.62 + 8.05 + 1.46)"

      # --- PROCESS ENERGIES ---
      # Remove closing parenthesis
      energy_str <- gsub(")", "", energy_part, fixed = TRUE)

      # Split by " = "
      # Left: Binding Energy, Right: Components
      e_split <- strsplit(energy_str, " = ", fixed = TRUE)[[1]]
      dG_binding <- as.numeric(e_split[1])

      # Split Components by " + "
      comps <- as.numeric(strsplit(e_split[2], " + ", fixed = TRUE)[[1]])

      # --- PROCESS COORDINATES ---
      # Split coords_part by " : " to separate Target vs Guide
      c_split <- strsplit(coords_part, " : ", fixed = TRUE)[[1]]
      target_side <- c_split[1]
      guide_side  <- c_split[2]

      # Target Coords: Split by space, take the last element
      t_tokens <- strsplit(target_side, " ", fixed = TRUE)[[1]]
      t_coord_str <- t_tokens[length(t_tokens)] # e.g. "23,50"

      # Guide Coords: Just trim whitespace
      g_coord_str <- trimws(guide_side) # e.g. "6,33"

      # Split comma values
      t_vals <- as.numeric(strsplit(t_coord_str, ",", fixed = TRUE)[[1]])
      g_vals <- as.numeric(strsplit(g_coord_str, ",", fixed = TRUE)[[1]])

      results_list[[i]] <- data.frame(
        dG_binding     = dG_binding,
        dG_duplex      = comps[1],
        dG_open_target = comps[2],
        dG_open_spacer  = comps[3],
        Target_Start   = t_vals[1],
        Target_End     = t_vals[2],
        Guide_Start    = g_vals[1],
        Guide_End      = g_vals[2]
      )

    } else {
      # Handle cases where RNAup finds no significant interaction
      results_list[[i]] <- data.frame(
        dG_binding=NA, dG_duplex=NA, dG_open_target=NA, dG_open_spacer=NA,
        Target_Start=NA, Target_End=NA, Guide_Start=NA, Guide_End=NA
      )
    }
  }
  return(do.call(rbind, results_list))
}


#' Check Cas13a Guide Self-Folding
#'
#' @description Predicts the secondary structure of the full guide RNA
#' (Direct Repeat + Spacer). High internal structure in the spacer region
#' suggests poor performance.
#'
#' @param spacer_with_linker_seq A DNAStringSet or character vector of spacer sequences.
#' @param rnafold_path Path to RNAfold binary.
#'
#' @return A data.frame with MFE (Minimum Free Energy) and structure notation.
#' @keywords internal
#'
predict_spacer_self_structure <- function(spacer_with_linker_seq,
                                         rnafold_path = "RNAfold") {

  # Input for RNAfold
  # >Name
  # Sequence
  input_lines <- paste0(">", seq_along(spacer_with_linker_seq), "\n",
                        spacer_with_linker_seq, collapse = "\n")

  # Run RNAfold --noPS (no postscript files)
  args <- c("--noPS")
  out <- tryCatch({
    system2(command = rnafold_path, args = args,
            input = input_lines, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    out <- ""
  })

  # Output format:
  # >Name
  # Sequence
  # .((...)).  (-10.50)

  # We only care about the structure lines (every 3rd line starting from line 3)
  structure_lines <- out[seq(3, length(out), by = 3)]
  results_list <- lapply(structure_lines, function(line) {
    parts <- strsplit(line, " (", fixed = TRUE)[[1]]
    struct <- paste(parts[1:(length(parts)-1)], collapse = " (")
    energy_part <- parts[length(parts)]
    energy <- as.numeric(gsub(")", "", energy_part, fixed = TRUE))
    return(data.frame(Structure = struct, MFE = energy))
  })

  return(do.call(rbind, results_list))
}

#' Check 0MM off-targets for each guide
#' @description Check 0MM off-targets for each guide inside the list of
#' offtargets_seq.
#' @param guides_seq DNAStringSet of the guides
#' @param offtargets_seq list of DNAStringSet elements where we will look for
#' 0MM off-targets and sum up results
#' @return A vector with one number for each guide
#' @keywords internal
#'
get_offtarget_count_MM0 <- function(guides_seq, offtargets_seq) {
  ot_count <- rep(0, length(guides_seq))
  if (length(offtargets_seq) == 0) {
    return(ot_count)
  }
  for (i in seq_along(guides_seq)) {
    guide <- guides_seq[[i]]
    ot_count[i] <- sum(sapply(offtargets_seq, function(transcriptome) {
      mindex <- Biostrings::vmatchPattern(guide, transcriptome)
      nmatch_per_seq <- elementNROWS(mindex)
      sum(nmatch_per_seq)
    }))
  }
  return(ot_count)
}

# annotation = "../../genomes/Homo_sapiens/Homo_sapiens.sqlite"
# genome = BSgenome.Hsapiens.ENSEMBL.GRChg38.p14primaryassembly::BSgenome.Hsapiens.ENSEMBL.GRChg38.p14primaryassembly
#
# txdb <- if (tools::file_ext(annotation) %in% c("gff3", "gtf")) {
#   suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation))
# } else {
#   suppressWarnings(AnnotationDbi::loadDb(annotation))
# }

# library(Biostrings)
# library(GenomicFeatures)
# library(BSgenome)
# library(IRanges)
# offtargets_rds = c("../../genomes/Homo_sapiens/Homo_sapiens.transcriptome.rds")
# offtargets_fasta = c() # Empty for testing purpose but users will upload
# input_fasta = testthat::test_path("testdata", "kras.fa")
# position = 600
# primer3_path = "/home/ai/Soft/primer3/src/primer3_core"

#' Design twisted Cas13a Guide Self-Folding
#'
#' @description Input: Patient RNA (SARS-CoV-2) or Genomic DNA (KRAS).
#' Conversion/Amplification:
#'   If RNA: It is Reverse Transcribed into cDNA.
#' PCR Step: Primers amplify this DNA/cDNA.
#' Transcription (The Key Step): The PCR product is used as a template for T7 RNA Polymerase to create new RNA.
#' Detection: Cas13 detects this newly transcribed RNA.
#' Input defaults for tests
#'
#' @param input_fasta A DNAStringSet fasta file path.
#' @param position Position on input target to design for.
#' @param guide_length Normally 28
#' @param allowed_positions 3:5, where we want guide to have mutation
#' @param t7_promoter "GAAATTAATACGACTCACTATAGGG" we attach this to left primer
#' @param linker "GCGCT" linker sequence - this is default optimized 5-mer
#' @param primer3_path Path to primer3_core.
#' @param primers_per_guide 2, by default we design two primer pairs per guide
#' @param rnaup_path Path to RNAup binary.
#' @param rnafold_path Path to RNAfold binary.
#'
#' @return A data.frame with MFE (Minimum Free Energy) and structure notation.
#' @export
#'
design_twisted_cas13a <- function(
    input_fasta,
    position,
    guide_length = 28,
    allowed_positions = 3:5, # This means the target 'position' can fall on index 3:5 of the guide.
    t7_promoter = "GAAATTAATACGACTCACTATAGGG", # this is 5'-3'
    linker = "GCGCT", # this is 5'-3' and DNA!
    primer3_path = Sys.which("primer3_core"),
    primers_per_guide = 2,
    rnaup_path = Sys.which("RNAup"),
    rnafold_path = Sys.which("RNAfold")) {

  # This has T, U not allowed! - we will output U if required...
  target <- Biostrings::readDNAStringSet(input_fasta)[[1]]
  target_len <- length(target)
  # Calculate theoretical start coordinates for the guides
  potential_starts <- position - allowed_positions + 1
  potential_ends <- potential_starts + guide_length - 1
  # Filter for valid bounds
  # We must ensure the guide fits entirely within the transcript (1 to rna_len)
  # We reserve space for primers (18bp padding)
  min_bound <- 1 + 18
  max_bound <- target_len - 18
  valid_idx <- which(potential_starts >= min_bound & potential_ends <= max_bound)
  # Check if any guides are possible
  if (length(valid_idx) == 0) {
    stop("No valid guides and primers can be generated at this position due to transcript boundary constraints.")
  }

  # Generate the valid Ranges
  # We only select the starts that survived the boundary check
  valid_starts <- potential_starts[valid_idx]
  guide_target <- IRanges(start = valid_starts, width = guide_length)
  guide_target_views <- Views(target, guide_target)
  guide_target_seq <- methods::as(guide_target_views, "DNAStringSet")

  # Now we generate up to N primer pairs with T7 around our position
  primers <- design_primers_rna(
    target, position, guide_target_seq, primers_per_guide, primer3_path)
  if (nrow(primers) == 0) stop("Primer3 failed to find valid primers.")
  target_match <- match(primers$SEQUENCE_INTERNAL_OLIGO, as.character(guide_target_seq))
  # I think for Twisted Cas13a it on 5' end of target
  # in normal Cas13a its on 3' end of target
  primers$MISMATCH_POSITION <- allowed_positions[valid_idx][target_match]
  primers$SPACER_START <- start(guide_target)[target_match]
  primers$SPACER_LENGTH <- width(guide_target)[target_match]
  primers$PFS <- as.character(
    methods::as(Views(target, flank(guide_target, width = 1)), "DNAStringSet"))
  # Linker is 5'-3'
  linker <- Biostrings::RNAStringSet(DNAString(linker))
  # spacer is 5'-3'
  spacer_seq <- Biostrings::RNAStringSet(Biostrings::reverseComplement(DNAStringSet(primers$SEQUENCE_INTERNAL_OLIGO)))

  # # Off-target analysis for each of the guides, 0 MM on the fasta
  # # TODO for the moment this hits false-positives because of conflict with
  # # target sequence
  # offtargets_rds <- lapply(offtargets_rds, readRDS)
  # offtargets_transcripts <- lapply(offtargets_fasta, readDNAStringSet)
  # primers$GUIDE_OFFTARGET_COUNT <-
  #   get_offtarget_count_MM0(guides_seq, offtargets_rds) +
  #   get_offtarget_count_MM0(guides_seq, offtargets_transcripts)

  # Generate Full Amplicons & Dual Components
  primers$amp_seq_dna <- as.character(Biostrings::subseq(
    DNAStringSet(lapply(1:nrow(primers), function(x) target)),
    start = primers$PRIMER_FORWARD_START, width = primers$PRODUCT_SIZE))
  primers$PRIMER_FORWARD_SEQUENCE_WITH_T7 <- paste0(
    toupper(t7_promoter), toupper(primers$PRIMER_FORWARD_SEQUENCE))
  primers$SPACER_WITH_LINKER <- paste0(
    as.character(spacer_seq), as.character(linker))

  # Guide vs Amplicon (Simulated RNA transcript)
  # We check the binding of the spacer to the RNA produced by the amplicon
  # We need to run this for every ROW of primers
  rnaup <- predict_spacer_amplicon_interaction(
    spacer_with_linker_seq = primers$SPACER_WITH_LINKER,
    amp_seq_rna = RNAStringSet(DNAStringSet(primers$amp_seq_dna)),
    rnaup_path = rnaup_path,
    window_size = 30)


  # Run RNAfold on the dcrRNA (Variable part)
  rnafold <- predict_spacer_self_structure(
    spacer_with_linker_seq = primers$SPACER_WITH_LINKER,
    rnafold_path = rnafold_path)

  final_dt <- cbind(rnaup, rnafold, primers)
  final_dt$linker_is_open <-
    !stringr::str_detect(
      substr(final_dt$Structure,
             nchar(final_dt$Structure) - nchar(linker) + 1,
             nchar(final_dt$Structure)), "[()]")
  final_dt$MFE_is_over_minus8 <- final_dt$MFE > -8 # 4bp hairpin is allowed, but 5bp not


  # WHAT do we want to return and in which order?!
  # In thermodynamics lower values is stringer interaction.
  # G at 28 is decreasing efficiency - add a column like that - hard filter

  # Interpreting the RNAfold:
  # If we have () in the Structure - potential linker issue - hard filter
  # If we have very LOW MFE - closer to 0 is better (less stable self-structure)

  # Interpretung RNAup:
  # dG_binding - if we compare different numbers - (with mut vs no mut)
  # perfect matches should have LOWER value, the bigger the difference between
  # the two the stronger they are in discriminating variants
  # dG_duplex - The energy purely from the 28-nucleotide base-pairing, assuming the strands were already unfolded.
  # LOWER is better - as it gives stronger bind.
  # dG_open_spacer - Lower is better (less energy to open spacer)
  # dG_open_target - Lower is better (less energy to open target)
  final_dt <- final_dt[order(
    final_dt$PFS == "G",
    -final_dt$linker_is_open,
    -final_dt$MFE_is_over_minus8,
    final_dt$dG_binding,
    final_dt$dG_open_target
  ), ]
  sel_cols <- c("PFS", "linker_is_open", "MFE", "dG_binding", "dG_duplex", "dG_open_spacer", "dG_open_target",
                "MISMATCH_POSITION", "SPACER_WITH_LINKER", "SPACER_START", "SPACER_LENGTH","PRODUCT_SIZE", "PRIMER_FORWARD_SEQUENCE_WITH_T7", "PRIMER_REVERSE_SEQUENCE",
                "PRIMER_FORWARD_START", "PRIMER_FORWARD_LENGTH", "PRIMER_REVERSE_START", "PRIMER_REVERSE_LENGTH",
                "PRIMER_FORWARD_TM", "PRIMER_REVERSE_TM", "PRIMER_FORWARD_GC", "PRIMER_REVERSE_GC")
  final_dt <- final_dt[, sel_cols]
  return(final_dt)
}
