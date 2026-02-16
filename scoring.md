# Scoring and Ordering Logic

This is an algorithmic specification of how guide scoring, SNV scoring, SNV selection, and template ordering work in the HDR design pipeline.

## 1. Inputs and Control Knobs

Primary knobs that change ranking behavior:
- `optimization_scheme` in `{balanced, disruption_first, safety_first}`
- `maximum_mutations_per_template` (max introduced SNVs per template)
- `score_efficiency` (whether guide efficiency models are run)
- `alphagenome_context` (optional tissue/context filter for AlphaGenome)
- optional annotations: dbSNP, CADD, ClinVar, AlphaGenome

## 2. High-Level Pipeline

1. Normalize input variants.
2. Discover candidate guides on a variant-centered sequence.
3. Optionally score guides with multiple efficiency models and aggregate ranks.
4. Build candidate SNV positions from guide footprints on template sequence.
5. Generate all possible SNVs at those positions and filter invalid/risky loci.
6. Compute SNV scoring features and assign scheme-specific `priority_group`.
7. For each guide and each requested SNV count `mpt`:
   - sort SNVs by scheme
   - greedily select non-overlapping SNVs
   - require exact cardinality (`selected_count == mpt`)
8. Build template sequence and compute template metrics.
9. Compute `disruption_bin` and apply final template sort based on scheme.

## 3. Guide Discovery and Scoring

### 3.1 Guide discovery

- Cas9 logic:
  - forward PAM: `NGG`
  - reverse PAM: `CCN`
- Each guide record contains:
  - `original`: 20mer protospacer
  - `with_pam`: 23mer protospacer+PAM

### 3.2 Optional guide efficiency scoring

If `score_efficiency = TRUE`, run multiple scoring models per guide.
Current model set:
- Doench 2014
- Moreno-Mateos 2015
- Labuhn 2018

Scoring failures are tolerated per model per guide:
- failed model call => `NA` for that model/guide
- other models still contribute

### 3.3 Aggregated guide rank

For each successful model:
1. rank guides by descending model score.
2. for each guide, collect available ranks (ignore `NA`).
3. compute geometric mean of collected ranks.
4. `rank_by_scores` is rank of geometric mean (lower is better).

If no model succeeds for any guide, `rank_by_scores = NA` for all guides.

## 4. Candidate SNV Generation

Candidate SNV positions are derived from guide-template alignment, then constrained:

1. Map guide footprint positions on template.
2. Exclude homology-arm positions (only active region is eligible).
3. Convert footprint position to `position_in_guide` (strand-aware).
4. Remove positions that overlap original protected template variants.

At each surviving position, generate all 3 possible single-nucleotide substitutions.

## 5. Candidate SNV Filtering and Annotation

Apply filters in sequence:

1. Remove positions near splice-site windows.
2. Remove positions near CDS start/stop boundaries.
3. Optionally remove ClinVar-overlapping positions.
4. Evaluate coding consequence across affected transcripts.
   - Keep SNVs only if every affected transcript is either:
     - synonymous, or
     - moot due to already-broken transcript state, or
     - non-coding for that SNV.

Optional annotations attached to surviving SNVs:
- dbSNP matching information
- noncoding overlap flag
- CADD score
- AlphaGenome splice-effect features

## 6. SNV Feature Engineering

### 6.1 Disruption tier (`disruption_tier`, lower is better)

Derived from `position_in_guide`:
- tier 1: `position_in_guide >= 22`
- tier 2: `>= 17`
- tier 3: `>= 13`
- tier 4: `>= 10`
- tier 5: otherwise

### 6.2 dbSNP priority (`dbSNP_priority`, lower is better)

- 1: known compatible SNP allele exists
- 2: overlapping SNP exists but allele is not compatible
- 3: no informative SNP evidence

### 6.3 CADD feature

- `cadd_imputed = CADD`
- if missing: impute to `15`

### 6.4 AlphaGenome feature

Composite score:

`ag_composite = SPLICE_SITES + SPLICE_SITE_USAGE + (SPLICE_JUNCTIONS / 5)`

Processing behavior:
- use absolute effect magnitude
- optional context filtering by `alphagenome_context` before max aggregation
- missing values default to NA

### 6.5 Safety tier (`safety_tier`, lower is better)

Default is tier 3, then apply rules:

1. Tier 5 if noncoding overlap exists.
2. Tier 4 if predicted risky:
   - `cadd_imputed > benign_cadd_threshold` OR
   - `ag_composite > ag_threshold`
3. Tier 2 if real low-risk evidence:
   - real CADD present, low CADD, not risky, no noncoding overlap.
4. Tier 1 if known benign dbSNP and no noncoding overlap.

Default thresholds:
- `benign_cadd_threshold = 15`
- `ag_threshold = 1.5`

## 7. Scheme-Specific SNV Priority Group

`priority_group` is assigned from (`disruption_tier`, `safety_tier`) and scheme.
Lower group is better.

### 7.1 `balanced`

- group 1: `d<=2 & s<=2`
- group 2: `d=3 & s<=2`
- group 3: `d<=2 & s=3`
- group 4: `d=3 & s=3`
- group 5: `d>=4 & s<=2`
- group 6: `d>=4 & s=3`
- group 7: `d<=2 & s=4`
- group 8: `s>=4 & d>=3`
- group 9: `s=5`

### 7.2 `disruption_first`

- group 1: `d<=2 & s<=4`
- group 2: `d=3 & s<=4`
- group 3: `d>=4 & s<=4`
- group 4: `s=5`

### 7.3 `safety_first`

- group 1..5 maps directly from `safety_tier` 1..5.

## 8. Per-Guide SNV Sorting and Selection

For each guide and each requested SNV count `mpt`:

### 8.1 Sort candidate SNVs by scheme

Ascending order unless noted:

- `balanced`:
  - `priority_group`
  - `disruption_tier`
  - `cadd_imputed`
  - `ag_impact_score`
  - `-position_in_guide` (higher position first)

- `disruption_first`:
  - `priority_group`
  - `safety_tier`
  - `cadd_imputed`
  - `ag_impact_score`

- `safety_first`:
  - `priority_group`
  - `disruption_tier`
  - `ag_impact_score`
  - `-position_in_guide`

### 8.2 Greedy non-overlap selection

Traverse sorted SNVs and accept an SNV only if both are true:
- no genomic overlap with selected SNVs
- no codon-key collision with selected SNVs (`tx_id + codon_num`)

Stop at `N = mpt` or exhaustion.

### 8.3 Exact-cardinality requirement

If selected SNVs count is not exactly `mpt`, that `(guide, mpt)` template is skipped.

## 9. Template Metrics

For each produced template:

- `snvs_introduced`: semicolon list of selected SNV ids
- `total_cadd`: sum of selected `cadd_imputed`
- `total_snp_quality_score`: sum of selected `dbSNP_priority`
- `max_alphagenome_score`: max selected AG composite
- `any_overlaps_noncoding`: any selected SNV overlaps noncoding annotations
- disruption metrics from guide-vs-template alignment:
  - `pam_disrupted_count`
  - `seed_disrupted_count`
  - `total_disruption_count`
- alignment strings for inspection:
  - `aln_guide`
  - `aln_template`

Seed/PAM mismatch windows are strand-aware:
- plus strand:
  - PAM indices: 22-23
  - seed indices: 11-20
- minus strand:
  - PAM indices: 1-2
  - seed indices: 4-13

## 10. Disruption Bin and Final Template Ordering

### 10.1 `disruption_bin` (lower is better)

Base bins use PAM status and total disruption count:

- bin 0 (Platinum): `pam_disrupted_count >= 1` AND `total_disruption_count >= 2`
- bin 1 (Gold): `pam_disrupted_count >= 1` AND `total_disruption_count == 1`
- bin 2 (Silver): `pam_disrupted_count == 0` AND `total_disruption_count >= 3`
- bin 3 (Bronze): `pam_disrupted_count == 0` AND `total_disruption_count == 2`
- bin 4 (Iron): `pam_disrupted_count == 0` AND `total_disruption_count == 1`
- bin 5 (Dirt): `total_disruption_count == 0`

Safety penalty is then applied directly to disruption rank
(except in `disruption_first`, where penalty is not applied):

- `unsafe_snv_count = number of selected SNVs with safety_tier >= 4`
- final `disruption_bin = min(5, base_disruption_bin + unsafe_snv_count)`

`n_snvs` is derived from `snvs_introduced` token count.
If `snvs_introduced` is empty, `n_snvs = 0`.

### 10.2 Final ordering per `optimization_scheme`

- `balanced`:
  1. `any_overlaps_noncoding` (FALSE first)
  2. `disruption_bin`
  3. `total_snp_quality_score`
  4. `n_snvs`
  5. `max_alphagenome_score`

- `disruption_first`:
  1. `disruption_bin`
  2. `n_snvs`
  3. `total_disruption_count`

- `safety_first`:
  1. `any_overlaps_noncoding`
  2. `total_snp_quality_score`
  3. `disruption_bin`
  4. `n_snvs`

## 11. Determinism and Edge Behavior

Deterministic behavior relies on:
- normalized/sorted variants
- deterministic sort keys at SNV and template stages
- fixed random seed for any random components in the broader workflow

Important edge behavior:
- no valid SNV positions in active region => stop
- no valid candidate SNVs after filtering => stop
- no successful guide scorers => guide rank remains `NA`
- missing CADD => imputed neutral value
- missing AlphaGenome modalities => treated as zero contribution
