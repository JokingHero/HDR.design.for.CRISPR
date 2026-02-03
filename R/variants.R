#' Check for VCF compliance (Pure Predicate)
#'
#' Returns TRUE if the variant is INVALID.
#' Does NOT throw errors (so we can use it to identify rows to fix).
#'
#' @param variants GRanges object
#' @return Logical vector (TRUE = Invalid/Needs Fixing, FALSE = Good)
#'
is_vcf_invalid <- function(variants) {
  if (length(variants) == 0) return(logical(0))

  ref <- as.character(variants$REF)
  alt <- as.character(variants$ALT)

  # 1. Check for Empty Alleles (e.g. REF="A", ALT="")
  has_empty_allele <- (nchar(ref) == 0) | (nchar(alt) == 0)

  # 2. Check for Unanchored Indels (e.g. REF="CT", ALT="A")
  # Condition: Lengths differ AND First bases do not match
  is_indel <- nchar(ref) != nchar(alt)

  # We initialize unanchored as FALSE
  unanchored_indel <- rep(FALSE, length(variants))

  # Only check variants that are Indels AND NOT empty (empty handled by #1)
  check_idx <- is_indel & !has_empty_allele

  if (any(check_idx)) {
    # Compare first character of REF vs ALT
    first_ref <- substr(ref[check_idx], 1, 1)
    first_alt <- substr(alt[check_idx], 1, 1)
    unanchored_indel[check_idx] <- (first_ref != first_alt)
  }

  return(has_empty_allele | unanchored_indel)
}

#' The Gatekeeper: Strict Assertion
#'
#' This runs AFTER correction. If this fails, the software stops.
#' @param variants Your variants to asses
#'
assert_vcf_valid <- function(variants) {
  invalid_mask <- is_vcf_invalid(variants)

  if (any(invalid_mask)) {
    bad_vars <- variants[invalid_mask]
    ref <- as.character(bad_vars$REF)
    alt <- as.character(bad_vars$ALT)
    is_empty <- (nchar(ref) == 0) | (nchar(alt) == 0)

    msg <- "VCF VALIDATION FAILED:"
    if (any(is_empty)) {
      msg <- paste(msg, "\nEmpty alleles detected (Recovery failed):",
                   head(paste(seqnames(bad_vars[is_empty]), start(bad_vars[is_empty]), sep=":")))
    }
    if (any(!is_empty)) {
      msg <- paste(msg, "\nUnanchored indels detected (First base mismatch):",
                   head(paste(seqnames(bad_vars[!is_empty]), start(bad_vars[!is_empty]), sep=":")))
    }

    stop(msg)
  }
}

#' Correct non-VCF compliant variants by adding genomic anchors.
#'
#' This function uses `is_vcf_invalid` to identify variants that are non-compliant
#' due to empty alleles (e.g., REF="A", ALT=""). It corrects them by prepending
#' the preceding genomic base to both REF and ALT and adjusting the coordinates.
#'
#' @param variants A GRanges object of variants to check and correct.
#' @param genome A BSgenome object or similar for use with getSeq.
#'
#' @return A GRanges object where all validatable variants are now VCF-compliant.
#' @import GenomicRanges
#' @import IRanges
#' @import Biostrings
#' @importFrom S4Vectors mcols
#'
correct_variants <- function(variants, genome) {
  invalid_mask <- is_vcf_invalid(variants)
  if (!any(invalid_mask)) {
    return(variants)
  }

  vars_to_fix <- variants[invalid_mask]
  anchor_gr <- flank(vars_to_fix, width = 1, start = TRUE)
  if (any(start(anchor_gr) < 1)) {
    # We only proceed with the ones we CAN fix
    fixable_idx <- start(anchor_gr) >= 1
    if (!any(fixable_idx)) return(variants)
    vars_to_fix <- vars_to_fix[fixable_idx]
    anchor_gr <- anchor_gr[fixable_idx]
    invalid_mask[invalid_mask] <- fixable_idx
  }
  anchor_bases <- as.character(getSeq(genome, anchor_gr))
  new_ref <- paste0(anchor_bases, as.character(vars_to_fix$REF))
  new_alt <- paste0(anchor_bases, as.character(vars_to_fix$ALT))
  ranges(vars_to_fix) <- IRanges(
    start = start(vars_to_fix) - 1,
    width = nchar(new_ref))

  vars_to_fix$REF <- new_ref
  vars_to_fix$ALT <- new_alt
  variants[invalid_mask] <- vars_to_fix
  return(variants)
}

#' Find variants that are not represented in a VCF-compliant manner.
#'
#' This function checks for two common VCF representation rules:
#' 1. REF and ALT alleles must not be empty strings.
#' 2. For indels (where nchar(REF) != nchar(ALT)), the first base of REF and ALT
#'    must be identical.
#' @param variants A GRanges object representing variants. Must have 'REF'
#'   and 'ALT' metadata columns.
#' @return An integer vector of indices corresponding to non-compliant variants.
#' @import GenomicRanges
#'
is_vcf_invalid <- function(variants) {
  if (!all(c("REF", "ALT") %in% names(mcols(variants)))) {
    stop("Input GRanges must have 'REF' and 'ALT' metadata columns.")
  }
  if (length(variants) == 0) {
    return(integer(0))
  }
  ref <- variants$REF
  alt <- variants$ALT
  has_empty_allele <- (nchar(ref) == 0) | (nchar(alt) == 0)
  is_indel <- nchar(ref) != nchar(alt)

  unanchored_indel <- rep(FALSE, length(variants))
  if (any(is_indel)) {
    maybe_una <- is_indel & !has_empty_allele
    unanchored_indel[maybe_una] <-
      substr(ref[maybe_una], 1, 1) != substr(alt[maybe_una], 1, 1)
  }
  return(has_empty_allele | unanchored_indel)
}

#' Validate and normalize variants against a reference genome.
#'
#' This function performs several checks and corrections:
#' 1. Validates that the REF allele in the input matches the reference genome.
#' 2. Detects and corrects for variants specified on the reverse strand.
#' 3. Normalizes non-VCF compliant indels (e.g., ALT="") by adding a
#'    padding base from the genome.
#' 4. Throws an error for unrecoverable issues like incorrect REF or
#'    ambiguous unanchored indels.
#'
#' @param variants_genomic A GRanges object of variants with REF and ALT columns.
#' @param genome A BSgenome object or similar for use with getSeq.
#'
#' @return A sorted, named, and normalized GRanges object ready for downstream use.
#' @import GenomicRanges
#' @import Biostrings
#' @import BSgenome
#'
normalize_variants <- function(variants_genomic, genome) {
  if (length(variants_genomic) == 0) return(variants_genomic)

  if (!all(c("REF", "ALT") %in% names(mcols(variants_genomic)))) {
    stop("Input GRanges must have 'REF' and 'ALT' metadata columns.")
  }
  if (any(width(variants_genomic) > 50)) {
    stop("We don't support variants larger than 50bp.")
  }
  genome_ref <- getSeq(genome, variants_genomic)
  user_ref <- DNAStringSet(variants_genomic$REF)

  # Check for a reverse-complement match first
  user_ref_rc <- reverseComplement(user_ref)
  if (all(genome_ref == user_ref_rc)) {
    message("NOTE: Provided REF alleles match the reverse-complement of the genome. Correcting REF and ALT to the plus strand.")
    variants_genomic$REF <- as.character(user_ref_rc)
    variants_genomic$ALT <- as.character(reverseComplement(DNAStringSet(variants_genomic$ALT)))
    user_ref <- DNAStringSet(variants_genomic$REF)
  }

  matches <- genome_ref == user_ref
  if (!all(matches)) {
    stop(c("REF allele mismatch detected for the following variants:",
           paste("\n", as.character(user_ref)[!matches],
                 "our genome has:",
                 as.character(genome_ref)[!matches])))
  }

  variants_genomic <- correct_variants(variants_genomic, genome)
  assert_vcf_valid(variants_genomic)
  variants_genomic <- sort(variants_genomic)
  names(variants_genomic) <- paste0("Variant ", seq_along(variants_genomic))
  return(variants_genomic)
}
