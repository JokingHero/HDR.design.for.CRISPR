#!/usr/bin/env python3
"""
CRISPR-MFH Scoring Script

Supports two modes:
  - offtarget: Batch-processes a CHOPOFF annotated CSV to score off-targets.
               Appends hardcoded PAM sequences and handles bulge trimming.
  - template:  Scores guide-vs-HDR-template pairs where PAM is already
               embedded in the alignment strings (from R pipeline).
               Handles strand orientation (reverse complement for - strand).

Usage:
    # Off-target mode (default):
    python3 predict.py --mode offtarget --input offtargets.csv --output scored.csv

    # Template mode (HDR pipeline integration):
    python3 predict.py --mode template --input templates.csv --output scored.csv

The model path is auto-resolved relative to this script's directory.
"""

import os
# Suppress TensorFlow logging before importing
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import sys
import csv
import argparse
import numpy as np
import tensorflow as tf


# ─── Constants ───────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_PATH = os.path.join(SCRIPT_DIR, 'model')
TLEN = 24  # Model expects 24-length sequences
GUIDE_PAM = 'NGG'  # Appended to guide alignment (offtarget mode only)
OFFTARGET_PAM = 'AGG'  # Appended to off-target alignment (offtarget mode only)
BATCH_SIZE = 10000  # TF inference batch size

COMPLEMENT = str.maketrans('ATGCatgc-', 'TACGtacg-')


# ─── Encoding ────────────────────────────────────────────────────────────────

CODE_DICT = {
    'A': [1, 0, 0, 0, 0],
    'T': [0, 1, 0, 0, 0],
    'G': [0, 0, 1, 0, 0],
    'C': [0, 0, 0, 1, 0],
    '-': [0, 0, 0, 0, 1],
}
DIRECTION_DICT = {'A': 5, 'G': 4, 'C': 3, 'T': 2, '-': 1}


def MFH_encoding(guide_seq, off_seq):
    """
    Encodes a single guide and off-target pair into the multi-feature matrices
    expected by the CRISPR-MFH model.
    """
    # Pad or truncate to exactly TLEN
    guide_seq = "-" * (TLEN - len(guide_seq)) + guide_seq.upper()
    off_seq = "-" * (TLEN - len(off_seq)) + off_seq.upper()

    if len(guide_seq) > TLEN:
        guide_seq = guide_seq[-TLEN:]
    if len(off_seq) > TLEN:
        off_seq = off_seq[-TLEN:]

    g_list = list(guide_seq)
    o_list = list(off_seq)

    pair_code = []
    on_mat = np.zeros((TLEN, 5), dtype=np.float32)
    off_mat = np.zeros((TLEN, 5), dtype=np.float32)

    for i in range(len(g_list)):
        # Resolve ambiguous bases
        if g_list[i] == 'N':
            g_list[i] = o_list[i]
        if g_list[i] == '_':
            g_list[i] = '-'
        if o_list[i] == '_':
            o_list[i] = '-'

        g_vec = np.array(CODE_DICT.get(g_list[i], [0, 0, 0, 0, 0]))
        o_vec = np.array(CODE_DICT.get(o_list[i], [0, 0, 0, 0, 0]))

        diff_code = np.maximum(g_vec, o_vec)

        dir_code = np.zeros(2)
        val_g = DIRECTION_DICT.get(g_list[i], 0)
        val_o = DIRECTION_DICT.get(o_list[i], 0)

        if g_list[i] != '-' and o_list[i] != '-' and val_g != val_o:
            if val_g > val_o:
                dir_code[0] = 1
            else:
                dir_code[1] = 1

        pair_code.append(np.concatenate((diff_code, dir_code)))
        on_mat[i] = g_vec
        off_mat[i] = o_vec

    return np.array(pair_code), on_mat, off_mat


# ─── Reverse Complement ─────────────────────────────────────────────────────

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence string.
    Handles A, T, G, C, and gap '-' characters.
    """
    return seq.translate(COMPLEMENT)[::-1]


# ─── Bulge Trimming ─────────────────────────────────────────────────────────

def trim_to_max_one_bulge(guide_aln, ref_aln):
    """
    Trim positions from the 3' end (right side) of both aligned sequences
    until at most 1 gap character ('-') remains across both sequences combined.

    CHOPOFF allows up to 3 bulges, but CRISPR-MFH only supports 1.
    This is called BEFORE PAM appending.
    """
    while (guide_aln.count('-') + ref_aln.count('-')) > 1:
        if len(guide_aln) <= 1:
            break
        guide_aln = guide_aln[:-1]
        ref_aln = ref_aln[:-1]
    return guide_aln, ref_aln


# ─── Preprocessing ───────────────────────────────────────────────────────────

def prepare_pair(guide_aln, ref_aln):
    """
    Prepare a guide/off-target alignment pair for MFH model input:
    1. Trim to max 1 bulge (if needed)
    2. Append PAM sequences
    """
    total_gaps = guide_aln.count('-') + ref_aln.count('-')
    if total_gaps > 1:
        guide_aln, ref_aln = trim_to_max_one_bulge(guide_aln, ref_aln)

    # Append PAM
    guide_with_pam = guide_aln + GUIDE_PAM
    ref_with_pam = ref_aln + OFFTARGET_PAM

    return guide_with_pam, ref_with_pam


# ─── Batch Prediction ───────────────────────────────────────────────────────

def predict_batch(model_infer, sequence_pairs):
    """
    Run MFH model inference on a batch of (guide, off-target) pairs.

    Args:
        model_infer: The TF model's serving function
        sequence_pairs: List of (guide_seq, off_target_seq) tuples (with PAM already appended)

    Returns:
        np.ndarray: 1D array of cleavage probability scores
    """
    if not sequence_pairs:
        return np.array([])

    batch_p, batch_on, batch_off = [], [], []
    for guide, off_target in sequence_pairs:
        p, on, off = MFH_encoding(guide, off_target)
        batch_p.append(p)
        batch_on.append(on)
        batch_off.append(off)

    num_samples = len(sequence_pairs)
    X_p = np.array(batch_p, dtype=np.float32).reshape(num_samples, 1, TLEN, 7)
    X_on = np.array(batch_on, dtype=np.float32).reshape(num_samples, 1, TLEN, 5)
    X_off = np.array(batch_off, dtype=np.float32).reshape(num_samples, 1, TLEN, 5)

    output_dict = model_infer(
        input_predict=tf.constant(X_p),
        input_on=tf.constant(X_on),
        input_off=tf.constant(X_off),
    )

    # prob_class_1 = probability of being a true off-target (higher = more dangerous)
    predictions = output_dict['main_output'].numpy()[:, 1]
    return predictions


# ─── Template Mode Processing ────────────────────────────────────────────────

def process_template_csv(input_path, output_path, model_path):
    """
    Score guide-vs-template pairs from the HDR pipeline.

    Input CSV must have columns: aln_guide, aln_template, strand
    - aln_guide/aln_template are 23bp aligned strings WITH PAM already included
    - strand is "+" or "-"
    - For "+" strand (NGG): PAM is at the right (positions 21-23) → use as-is
    - For "-" strand (CCN): PAM is at the left (positions 1-3) → reverse complement

    No PAM appending or bulge trimming is performed in this mode.
    """
    # Load model once
    print(f"Loading CRISPR-MFH model from: {model_path}", file=sys.stderr)
    model = tf.saved_model.load(model_path)
    infer = model.signatures["serving_default"]
    print("Model loaded successfully.", file=sys.stderr)

    # Read all rows
    rows = []
    with open(input_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        original_fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    if not rows:
        print("No data rows found in input CSV. Writing empty output.", file=sys.stderr)
        with open(output_path, 'w', newline='') as f:
            if original_fieldnames:
                writer = csv.DictWriter(f, fieldnames=list(original_fieldnames) + ['crispr_mfh_score'])
                writer.writeheader()
        return

    print(f"Processing {len(rows)} template rows...", file=sys.stderr)

    # Prepare all sequence pairs
    prepared_pairs = []
    for row in rows:
        guide_seq = row.get('aln_guide', '').strip()
        template_seq = row.get('aln_template', '').strip()
        strand = row.get('strand', '+').strip()

        if not guide_seq or not template_seq:
            prepared_pairs.append(None)
            continue

        try:
            # For minus strand (CCN): reverse complement both to put PAM at 3' end
            if strand == '-':
                guide_seq = reverse_complement(guide_seq)
                template_seq = reverse_complement(template_seq)

            # Clamp to TLEN (24bp) — MFH model expects exactly 24-length inputs.
            # R pipeline produces 23bp alignments (20bp protospacer + 3bp PAM),
            # which MFH_encoding pads to 24 with a leading gap '-'.
            # If sequences are longer (e.g. with indels), truncate from the 5' end.
            if len(guide_seq) > TLEN:
                guide_seq = guide_seq[-TLEN:]
            if len(template_seq) > TLEN:
                template_seq = template_seq[-TLEN:]

            prepared_pairs.append((guide_seq.upper(), template_seq.upper()))
        except Exception as e:
            print(f"Warning: Failed to prepare template pair: {e}", file=sys.stderr)
            prepared_pairs.append(None)

    # Run inference in batches
    scores = [None] * len(rows)
    valid_indices = [i for i, p in enumerate(prepared_pairs) if p is not None]
    valid_pairs = [prepared_pairs[i] for i in valid_indices]

    for batch_start in range(0, len(valid_pairs), BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, len(valid_pairs))
        batch = valid_pairs[batch_start:batch_end]
        batch_indices = valid_indices[batch_start:batch_end]

        try:
            batch_scores = predict_batch(infer, batch)
            for idx, score in zip(batch_indices, batch_scores):
                scores[idx] = float(score)
        except Exception as e:
            print(f"Warning: Batch prediction failed (rows {batch_start}-{batch_end}): {e}", file=sys.stderr)

        if batch_end % (BATCH_SIZE * 5) == 0 or batch_end == len(valid_pairs):
            print(f"  Scored {batch_end}/{len(valid_pairs)} valid pairs...", file=sys.stderr)

    # Write output CSV
    output_fieldnames = list(original_fieldnames) + ['crispr_mfh_score']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames)
        writer.writeheader()
        for row, score in zip(rows, scores):
            row['crispr_mfh_score'] = f"{score:.6f}" if score is not None else ""
            writer.writerow(row)

    scored_count = sum(1 for s in scores if s is not None)
    print(f"Done. Scored {scored_count}/{len(rows)} templates. Output: {output_path}", file=sys.stderr)


# ─── CLI Batch Processing (Off-target mode) ──────────────────────────────────

def process_csv(input_path, output_path, model_path):
    """
    Read annotated CHOPOFF CSV, score each off-target with CRISPR-MFH,
    and write the result with an appended crispr_mfh_score column.
    """
    # Load model once
    print(f"Loading CRISPR-MFH model from: {model_path}", file=sys.stderr)
    model = tf.saved_model.load(model_path)
    infer = model.signatures["serving_default"]
    print("Model loaded successfully.", file=sys.stderr)

    # Read all rows first to determine fieldnames
    rows = []
    with open(input_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        original_fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    if not rows:
        print("No data rows found in input CSV. Writing empty output.", file=sys.stderr)
        with open(output_path, 'w', newline='') as f:
            if original_fieldnames:
                writer = csv.DictWriter(f, fieldnames=list(original_fieldnames) + ['crispr_mfh_score'])
                writer.writeheader()
        return

    print(f"Processing {len(rows)} off-target rows...", file=sys.stderr)

    # Prepare all sequence pairs
    prepared_pairs = []
    for row in rows:
        guide_aln = row.get('alignment_guide', '')
        ref_aln = row.get('alignment_reference', '')

        if not guide_aln or not ref_aln:
            prepared_pairs.append(None)
            continue

        try:
            guide_pam, ref_pam = prepare_pair(guide_aln, ref_aln)
            prepared_pairs.append((guide_pam, ref_pam))
        except Exception as e:
            print(f"Warning: Failed to prepare pair: {e}", file=sys.stderr)
            prepared_pairs.append(None)

    # Run inference in batches
    scores = [None] * len(rows)
    valid_indices = [i for i, p in enumerate(prepared_pairs) if p is not None]
    valid_pairs = [prepared_pairs[i] for i in valid_indices]

    for batch_start in range(0, len(valid_pairs), BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, len(valid_pairs))
        batch = valid_pairs[batch_start:batch_end]
        batch_indices = valid_indices[batch_start:batch_end]

        try:
            batch_scores = predict_batch(infer, batch)
            for idx, score in zip(batch_indices, batch_scores):
                scores[idx] = float(score)
        except Exception as e:
            print(f"Warning: Batch prediction failed (rows {batch_start}-{batch_end}): {e}", file=sys.stderr)
            # Leave scores as None for this batch

        if batch_end % (BATCH_SIZE * 5) == 0 or batch_end == len(valid_pairs):
            print(f"  Scored {batch_end}/{len(valid_pairs)} valid pairs...", file=sys.stderr)

    # Write output CSV
    output_fieldnames = list(original_fieldnames) + ['crispr_mfh_score']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames)
        writer.writeheader()
        for row, score in zip(rows, scores):
            row['crispr_mfh_score'] = f"{score:.6f}" if score is not None else ""
            writer.writerow(row)

    scored_count = sum(1 for s in scores if s is not None)
    print(f"Done. Scored {scored_count}/{len(rows)} off-targets. Output: {output_path}", file=sys.stderr)


# ─── Main ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CRISPR-MFH Scoring (Off-target and Template modes)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Off-target mode (CHOPOFF pipeline):
  python3 predict.py --mode offtarget --input offtargets.csv --output scored.csv

  # Template mode (HDR pipeline - guide vs template):
  python3 predict.py --mode template --input templates.csv --output scored.csv

  # Single pair test:
  python3 predict.py --test
        """,
    )
    parser.add_argument('--mode', choices=['offtarget', 'template'], default='offtarget',
                        help='Scoring mode: offtarget (CHOPOFF) or template (HDR pipeline)')
    parser.add_argument('--input', help='Input CSV file')
    parser.add_argument('--output', help='Output CSV file (with crispr_mfh_score column appended)')
    parser.add_argument('--model', default=DEFAULT_MODEL_PATH,
                        help='Path to TF SavedModel directory')
    parser.add_argument('--test', action='store_true',
                        help='Run a quick test with example sequences')

    args = parser.parse_args()

    if args.test:
        # Quick self-test
        print("--- CRISPR-MFH Self-Test ---")
        model = tf.saved_model.load(args.model)
        infer = model.signatures["serving_default"]

        test_pairs = [
            # 20bp guide + NGG PAM, 20bp off-target + AGG PAM
            ("GAGTCCGAGCAGAAGAAGAANGG", "GAGTCCGAGCAGAAGAAGAAGAG"),
        ]
        scores = predict_batch(infer, test_pairs)
        for pair, s in zip(test_pairs, scores):
            print(f"Guide:  {pair[0]}")
            print(f"OT:     {pair[1]}")
            print(f"Score:  {s:.5f}")
        print("--- Test complete ---")

    elif args.input and args.output:
        if args.mode == 'template':
            process_template_csv(args.input, args.output, args.model)
        else:
            process_csv(args.input, args.output, args.model)

    else:
        parser.print_help()
        sys.exit(1)