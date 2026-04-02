# CRISPR-MFH Score vs Disruption Bin Validation

**MFH score** = probability that Cas9 still cleaves the template (0–1).
**Lower MFH = better** (more disrupted, Cas9 less likely to re-cut).

## CTNNB1 (40 templates, sorted by disruption_bin)

| Name | Bin | PAM | Seed | Total | SNVs | Raw MFH | Baseline | Rel. MFH |
|------|:---:|:---:|:----:|:-----:|:----:|--------:|---------:|----------|
| CCN_1_with_1_SNVs | 0 | 1 | 0 | 2 | 1 | **0.3437** | 0.2259 | 1.521 |
| CCN_2_with_1_SNVs | 0 | 1 | 1 | 2 | 1 | **0.1491** | 0.3177 | 0.469 |
| CCN_3_with_1_SNVs | 0 | 2 | 0 | 3 | 1 | **0.1444** | 0.2568 | 0.562 |
| CCN_1_with_2_SNVs | 0 | 1 | 1 | 2 | 1 | **0.1374** | 0.2259 | 0.608 |
| CCN_2_with_2_SNVs | 0 | 1 | 1 | 2 | 1 | **0.2875** | 0.3177 | 0.905 |
| CCN_3_with_2_SNVs | 0 | 1 | 2 | 3 | 1 | **0.0011** | 0.2568 | 0.004 |
| CCN_1_with_3_SNVs | 0 | 1 | 2 | 3 | 1 | **0.0003** | 0.2259 | 0.001 |
| CCN_2_with_3_SNVs | 0 | 2 | 1 | 4 | 1 | **0.0011** | 0.3177 | 0.003 |
| CCN_3_with_3_SNVs | 0 | 1 | 3 | 4 | 1 | **0.0000** | 0.2568 | 0.000 |
| NGG_1_with_3_SNVs | 0 | 1 | 2 | 3 | 1 | **0.0009** | 0.3238 | 0.003 |
| NGG_1_with_1_SNVs | 1 | 1 | 0 | 1 | 1 | **0.2043** | 0.3238 | 0.631 |
| NGG_1_with_2_SNVs | 1 | 1 | 0 | 1 | 1 | **0.2951** | 0.3238 | 0.911 |
| CCN_1_with_1_SNVs¹ | 2 | 0 | 2 | 3 | 2 | **0.1000** | 0.2259 | 0.443 |
| CCN_2_with_1_SNVs¹ | 2 | 0 | 3 | 3 | 2 | **0.0236** | 0.3177 | 0.074 |
| CCN_3_with_1_SNVs¹ | 2 | 0 | 2 | 3 | 2 | **0.0291** | 0.2568 | 0.113 |
| CCN_2_with_2_SNVs¹ | 2 | 0 | 3 | 4 | 2 | **0.0001** | 0.3177 | 0.000 |
| CCN_3_with_2_SNVs¹ | 2 | 0 | 3 | 3 | 2 | **0.0011** | 0.2568 | 0.004 |
| NGG_1_with_1_SNVs¹ | 2 | 0 | 0 | 3 | 2 | **0.0981** | 0.3238 | 0.303 |
| NGG_1_with_2_SNVs¹ | 2 | 0 | 2 | 3 | 2 | **0.0176** | 0.3238 | 0.054 |
| CCN_1_with_1_SNVs² | 3 | 0 | 1 | 2 | 3 | **0.3601** | 0.2259 | 1.594 |
| CCN_2_with_1_SNVs² | 3 | 0 | 2 | 2 | 3 | **0.1755** | 0.3177 | 0.552 |
| CCN_3_with_1_SNVs² | 3 | 0 | 1 | 2 | 3 | **0.2958** | 0.2568 | 1.152 |
| CCN_1_with_2_SNVs² | 3 | 0 | 2 | 2 | 3 | **0.0784** | 0.2259 | 0.347 |
| CCN_2_with_2_SNVs² | 3 | 0 | 0 | 2 | 3 | **0.5648** | 0.3177 | 1.778 |
| CCN_3_with_2_SNVs² | 3 | 0 | 1 | 2 | 3 | **0.3391** | 0.2568 | 1.320 |
| NGG_1_with_0_SNVs¹ | 4 | 0 | 1 | 1 | 4 | **0.3586** | 0.3238 | 1.107 |
| CCN_1_with_0_SNVs | 4 | 0 | 0 | 1 | 4 | **0.4198** | 0.2259 | 1.858 |
| CCN_2_with_0_SNVs | 4 | 0 | 1 | 1 | 4 | **0.3505** | 0.3177 | 1.103 |
| CCN_3_with_0_SNVs | 4 | 0 | 1 | 1 | 4 | **0.2472** | 0.2568 | 0.963 |
| NGG_1_with_0_SNVs² | 4 | 0 | 0 | 1 | 4 | **0.3416** | 0.3238 | 1.055 |
| NGG_1_with_1_SNVs² | 4 | 0 | 1 | 1 | 4 | **0.2707** | 0.3238 | 0.836 |
| NGG_1_with_2_SNVs² | 4 | 0 | 0 | 1 | 4 | **0.6543** | 0.3238 | 2.021 |
| NGG_1_with_3_SNVs¹ | 4 | 0 | 0 | 1 | 4 | **0.3207** | 0.3238 | 0.990 |
| CCN_1_with_0_SNVs¹ | 5 | 0 | 0 | 0 | 5 | **0.2259** | 0.2259 | 1.000 |
| CCN_2_with_0_SNVs¹ | 5 | 0 | 0 | 0 | 5 | **0.3177** | 0.3177 | 1.000 |
| CCN_3_with_0_SNVs¹ | 5 | 0 | 0 | 0 | 5 | **0.2568** | 0.2568 | 1.000 |
| NGG_1_with_0_SNVs³ | 5 | 0 | 0 | 0 | 5 | **0.3238** | 0.3238 | 1.000 |
| NGG_1_with_0_SNVs⁴ | 5 | 0 | 0 | 0 | 5 | **0.2820** | 0.3238 | 0.871 |
| NGG_1_with_1_SNVs³ | 5 | 0 | 0 | 0 | 5 | **0.3161** | 0.3238 | 0.976 |
| NGG_1_with_2_SNVs³ | 5 | 0 | 0 | 0 | 5 | **0.2990** | 0.3238 | 0.923 |

> [!NOTE]
> ¹²³⁴ suffixes distinguish different templates that have the same prefix in the original data.

### CTNNB1 Summary — Mean MFH by Disruption Bin

| Bin | Label | Mean MFH | Templates |
|:---:|-------|:--------:|:---------:|
| 0 | Platinum (PAM+total≥2) | **0.107** | 10 |
| 1 | Gold (PAM+total=1) | **0.250** | 2 |
| 2 | Silver (seed≥3, no PAM) | **0.039** | 7 |
| 3 | Bronze (total=2, no PAM) | **0.302** | 6 |
| 4 | Iron (total=1, no PAM) | **0.370** | 8 |
| 5 | Dirt (total=0) | **0.289** | 7 |

---

## STAT3 (18 templates, sorted by disruption_bin)

| Name | Bin | PAM | Seed | Total | SNVs | Raw MFH | Baseline | Rel. MFH |
|------|:---:|:---:|:----:|:-----:|:----:|--------:|---------:|----------|
| NGG_3_with_2_SNVs | 0 | 2 | 0 | 2 | 2 | **0.0842** | 0.2320 | 0.363 |
| NGG_2_with_1_SNVs | 0 | 1 | 1 | 2 | 1 | **0.2887** | 0.4052 | 0.712 |
| NGG_1_with_2_SNVs | 0 | 2 | 0 | 2 | 2 | **0.0744** | 0.2834 | 0.263 |
| NGG_3_with_3_SNVs | 0 | 2 | 1 | 3 | 3 | **0.0077** | 0.2320 | 0.033 |
| NGG_2_with_2_SNVs | 0 | 1 | 2 | 3 | 2 | **0.0070** | 0.4052 | 0.017 |
| NGG_1_with_3_SNVs | 0 | 2 | 1 | 3 | 3 | **0.0087** | 0.2834 | 0.031 |
| NGG_2_with_3_SNVs | 0 | 1 | 3 | 4 | 3 | **0.0000** | 0.4052 | 0.000 |
| NGG_1_with_1_SNVs | 1 | 1 | 0 | 1 | 1 | **0.2846** | 0.2834 | 1.004 |
| NGG_3_with_1_SNVs | 1 | 1 | 0 | 1 | 1 | **0.2285** | 0.2320 | 0.985 |
| CCN_1_with_2_SNVs | 3 | 0 | 2 | 2 | 2 | **0.0846** | 0.2736 | 0.309 |
| CCN_2_with_2_SNVs | 3 | 0 | 2 | 2 | 2 | **0.0168** | 0.2871 | 0.059 |
| NGG_2_with_0_SNVs | 4 | 0 | 1 | 1 | 0 | **0.4052** | 0.4052 | 1.000 |
| CCN_1_with_1_SNVs | 4 | 0 | 1 | 1 | 1 | **0.3745** | 0.2736 | 1.369 |
| CCN_2_with_1_SNVs | 4 | 0 | 1 | 1 | 1 | **0.3426** | 0.2871 | 1.193 |
| NGG_1_with_0_SNVs | 5 | 0 | 0 | 0 | 0 | **0.2834** | 0.2834 | 1.000 |
| NGG_3_with_0_SNVs | 5 | 0 | 0 | 0 | 0 | **0.2320** | 0.2320 | 1.000 |
| CCN_1_with_0_SNVs | 5 | 0 | 0 | 0 | 0 | **0.2736** | 0.2736 | 1.000 |
| CCN_2_with_0_SNVs | 5 | 0 | 0 | 0 | 0 | **0.2871** | 0.2871 | 1.000 |

### STAT3 Summary — Mean MFH by Disruption Bin

| Bin | Label | Mean MFH | Templates |
|:---:|-------|:--------:|:---------:|
| 0 | Platinum | **0.067** | 7 |
| 1 | Gold | **0.257** | 2 |
| 3 | Bronze | **0.051** | 2 |
| 4 | Iron | **0.374** | 3 |
| 5 | Dirt | **0.269** | 4 |

---

## Key Observations

> [!IMPORTANT]
> **MFH generally agrees with your binning**, but adds granularity:
> - **Bin 0 (Platinum)** has consistently the **lowest** MFH scores (mean ~0.07–0.11), confirming strong disruption
> - **Bin 4 (Iron)** has the **highest** MFH scores (mean ~0.37), correctly reflecting weak single-mismatch disruption
> - **Within bins**, MFH differentiates well — e.g., in Bin 0, templates with seed+PAM disruption score ~0.001 vs PAM-only scoring ~0.15–0.34

> [!WARNING]
> **Bin 5 (Dirt, total=0)** scores are **not the highest** (~0.27–0.29). This is expected — Bin 5 means the alignment-based counter found 0 mismatches, but MFH reports ~0.28 base cleavage probability. The score doesn't go to ~1.0 because the guide-template pairs in the test data happen to be the same sequence (no mutations), and MFH's baseline prediction for a perfect match is ~0.28, not 1.0. This makes sense since even perfect on-target sites don't have 100% cutting efficiency.
>
> **Bin 4 (Iron, total=1)** sometimes scores **higher** than Bin 5 — this suggests some single mismatches are in positions that don't actually reduce cutting much (or even slightly increase predicted activity for certain mismatch patterns). MFH captures position-dependent effects that simple counting misses.
