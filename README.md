# HDR.design.for.CRISPR

**HDR design for CRISPR**

`HDR.design.for.CRISPR` is an R package designed to automate the creation of Homology-Directed Repair (HDR) templates for CRISPR experiments. It generates repair templates using various strategies (balanced, safety-first, disruption-first), designs qPCR probes for validation, and integrates state-of-the-art scoring methods including **DeepMind's AlphaGenome** and **CADD**.

---

## ⚠️ Prerequisites & System Dependencies

**Important:** Before installing this R package, you must ensure the following system dependencies and Python environments are set up.

### 1. Python Environment & AlphaGenome (Optional)
This package relies on Google DeepMind's AlphaGenome for variant scoring.

1.  **Get your API Key:**
    You must obtain your own access key from the official site:
    [https://deepmind.google.com/science/alphagenome/](https://deepmind.google.com/science/alphagenome/)

2.  **Install Python Dependencies:**
    Ensure you have **Python 3** installed. Then, run the following command in your terminal to install the required libraries:

    ```bash
    pip install -U alphagenome
    pip install -U argparse
    pip install -U pandas
    ```

### 2. Primer3 (Optional)  

To utilize the PCR primer design functionality, you need `primer3`.

* **MacOS (via Homebrew):** `brew install primer3`
* **Linux (via apt):** `sudo apt-get install primer3`
* **Manual:** Download the `primer3_core` executable from [primer3-org](https://github.com/primer3-org/primer3).

---

## 📦 Installation

### Step 1: Install Bioconductor Dependencies
This package depends on several genomic data packages. It is recommended to install these via `BiocManager` first to ensure all binary dependencies are resolved correctly.  

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg38",
    "GenomicRanges",
    "GenomicFeatures",
    "VariantAnnotation",
    "crisprScore",
    "GenomicScores",
    "rtracklayer",
    "Biostrings",
    "pwalign",
    "txdbmaker"
))
```

### Step 2: Install HDR.design.for.CRISPR
You can install the development version from GitHub:  

```r
# install.packages("devtools")
devtools::install_github("JokingHero/HDR.design.for.CRISPR")
```

## 🚀 Usage

The core function of this package is `design_hdr`. Below is a functional example using the `STAT3` gene dataset included with the package.

### Example: STAT3 Mutation Correction

This example demonstrates how to design a repair template to **correct** a mutation.  

* **Scenario:** The cell line has a specific mutation (`A`) on the genome.  
* **Goal:** We want to revert it to the Wild Type (`G`).  

```r
library(HDR.design.for.CRISPR)
library(BSgenome.Hsapiens.UCSC.hg38)

# 1. Locate the example annotation file included in the package
# (This file is a mini-GFF3 provided for testing purposes)
annotation_file <- system.file("extdata", "gencode.v42.annotation_mini.gff3", 
                               package = "HDR.design.for.CRISPR")

# 2. Setup output directory
output_dir <- tempfile(pattern = "STAT3_design_")
dir.create(output_dir)

# 3. Run the design
# We are targeting the STAT3 gene on chr17 (minus strand).
design_hdr(
  design_name = "STAT3_correction",
  chrom = "chr17",
  variant_start = 42346635,
  variant_end = 42346635,
  REF = "G",   # The Wild Type base
  ALT = "A",   # The Mutation currently present in the cell
  
  # Logic Vectors:
  # TRUE: The cell line currently HAS the mutation (A)
  ALT_on_genome = TRUE,  
  
  # FALSE: We want the repair template to contain the REF (G), not the ALT
  # (Set this to TRUE if you were INTRODUCING a specific mutation instead)
  ALT_on_templates = FALSE,

  output_dir = output_dir,
  annotation = annotation_file,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  
  # Optimization settings
  optimization_scheme = "balanced",
  maximum_mutations_per_template = 3,
  
  # Optional: Improve design safety with known SNPs and CADD scores
  # (Requires installed 'SNPlocs' and 'cadd' packages)
  # snps = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
  # cadd = GenomicScores::getGScores("cadd.v1.6.hg38")
)

print(paste("Design results written to:", output_dir))
```

### Key Parameters Explained

| Parameter | Description |
| :--- | :--- |
| `optimization_scheme` | Controls how synonymous SNPs are selected. <br>• **balanced**: Balances safety, PAM disruption, and SNP quality.<br>• **safety_first**: Avoids non-coding overlaps/splice sites.<br>• **disruption_first**: Prioritizes PAM/Guide disruption. |
| `ALT_on_genome` | `TRUE/FALSE`. Does the cell line you are editing currently possess the ALT mutation? (Select `FALSE` if targeting Wild Type). |
| `ALT_on_templates` | `TRUE/FALSE`. Do you want the repair template to contain the ALT mutation? (Select `TRUE` if you are introducing the specific mutation). |
| `alphagenome_context` | (Optional) Cell types for which to filter AlphaGenome scores. |
| `clinvar` / `cadd` | (Optional) Paths to ClinVar VCF or CADD scoring functions to avoid introducing pathogenic synonymous mutations. |

---

## 📂 Output

The package generates a directory containing your design results. The output format and interpretation of scores are described in detail in the AlphaGenome documentation:

👉 **[Understanding Output & Variant Scoring](https://www.alphagenomedocs.com/variant_scoring.html)**

Typical outputs include:  
* **Repair Templates:** Sequences with optimal synonymous mutations.  
* **Guides:** Scored CRISPR guides (using `crisprScore`).  
* **Probes:** qPCR probes for HDR, NHEJ, and Reference detection.  
* **Primers:** PCR primers flanking the edit region (if `primer3` is provided).  

---

## ⚖️ License

AGPL (>= 3)
