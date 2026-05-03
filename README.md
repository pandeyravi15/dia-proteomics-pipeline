# Brain DIA Proteomics Analysis Pipeline

[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%3E%3D3.16-green)](https://bioconductor.org/)

A reproducible R-based pipeline for end-to-end analysis of **Data-Independent Acquisition (DIA) proteomics data** from mouse or human brain tissue. Originally developed for Alzheimer's disease mouse models but generalized for any case-control DIA proteomics study with multiple genotypes, sexes, and timepoints.

---

## 📋 What This Pipeline Does

| Step | Analysis |
|------|----------|
| QC | Protein abundance matrix loading, boxplot QC, summarization to protein level |
| PCA | Probabilistic PCA (handles missing values) with metadata overlay |
| Differential Abundance | ANOVA + Tukey HSD post-hoc, stratified by sex and age |
| Visualization | Volcano plots, DEP barplots, UpSet overlap plots, expression boxplots |
| Functional Enrichment | KEGG pathway ORA and GSEA via `clusterProfiler` |

---

## 🗂️ Repository Structure

```
dia-proteomics-pipeline/
│
├── README.md                          # This file
├── LICENSE                            # MIT license
├── .gitignore                         # Excludes data, results, and large files
│
├── notebooks/
│   └── DIA_Proteomics_Pipeline.Rmd    # Main analysis notebook (edit config at top)
│
├── R/
│   └── Helper_Functions_Prot.R        # Custom helper functions
│
├── data/
│   └── example/
│       ├── example_abundance.csv      # Example protein abundance matrix (50 proteins, 10 samples)
│       └── example_traits.csv         # Example sample metadata
│
├── figures/                           # Example output figures
│   ├── example_PCA.png
│   ├── example_volcano.png
│   └── example_KEGG.png
│
├── results/                           # gitignored — output saved here locally
│
└── docs/
    └── input_format.md                # Detailed input file format guide
```

---

## ⚙️ Requirements

### R Version
R ≥ 4.2.0

### CRAN Packages
```r
install.packages(c(
  "tidyverse", "gt", "UpSetR", "EnhancedVolcano", "data.table"
))
```

### Bioconductor Packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "pcaMethods",        # Probabilistic PCA
  "clusterProfiler",   # Pathway enrichment
  "AnnotationDbi",     # Gene ID mapping
  "org.Mm.eg.db",      # Mouse gene annotations (swap for org.Hs.eg.db for human)
  "enrichplot"         # Enrichment visualization
))
```

---

## 🚀 Quick Start

### Step 1: Prepare Your Input Files

You need two CSV files (see `docs/input_format.md` for details):

**1. Protein Abundance Matrix** (`post_TAMPOR_cleanDat.csv` or similar)
- Rows = proteins, Columns = samples
- Values should be pre-normalized (e.g., median-normalized, TAMPOR-corrected)
- ProteinID column can be `Protein` or `Protein|Peptide` format

**2. Sample Metadata** (`post_TAMPOR_traits.csv` or similar)
- One row per sample
- Required columns: Sample ID, Genotype, Sex, Age

### Step 2: Configure the Pipeline

Open `notebooks/DIA_Proteomics_Pipeline.Rmd` and edit the **USER CONFIGURATION** block at the top:

```r
# Key settings to change:
data_dir         <- "../data/"
abundance_file   <- "your_abundance_file.csv"
traits_file      <- "your_metadata_file.csv"

control_genotype <- "Your_Control_Genotype"
case_genotypes   <- c("Case_1", "Case_2", "Case_3")
age_groups       <- c(4, 12)        # your timepoints
organism_db      <- "org.Mm.eg.db"  # change to org.Hs.eg.db for human
kegg_organism    <- "mmu"           # "hsa" for human
```

### Step 3: Run the Analysis

In RStudio, open the notebook and click **Knit** to generate an HTML report, or run chunk by chunk interactively.

```r
rmarkdown::render("notebooks/DIA_Proteomics_Pipeline.Rmd")
```

---

## 📁 Input File Format

### Abundance Matrix

```
ProteinID,Sample1,Sample2,Sample3,...
MAPT|peptide1,12.3,11.8,12.1,...
APP|peptide2,10.5,10.2,10.8,...
```

- First column: ProteinID (can be `Protein` or `Protein|Peptide`)
- Remaining columns: one per sample, matching IDs in metadata

### Metadata File

```
X,Sex,Genotype,Age
Sample1,Female,Control,4
Sample2,Male,Control,4
Sample3,Female,CaseA,4
```

- `X`: Sample ID (must match column names in abundance matrix)
- `Sex`: Male / Female
- `Genotype`: group labels for each sample
- `Age`: numeric timepoint

See `docs/input_format.md` for more examples and edge cases.

---

## 🧪 Test with Example Data

```r
# Run with built-in example files
# In the USER CONFIGURATION block, set:
data_dir      <- "../data/example/"
abundance_file <- "example_abundance.csv"
traits_file    <- "example_traits.csv"
```

---

## 📊 Output

| Output | Location | Description |
|--------|----------|-------------|
| HTML Report | `notebooks/` | Full rendered analysis |
| ANOVA Results | `results/ANOVA_Results_DIA_Proteomics.Rdata` | All DEP results |
| Processed Data | `data/Processed_DIA_Proteomics.RData` | Cleaned abundance + traits |

---

## 🔬 Methods Summary

**Quality Control:** Sample distributions are assessed by boxplot after normalization. Proteins with pipe-delimited IDs are summarized to the protein level by averaging duplicate entries.

**Dimensionality Reduction:** Probabilistic PCA (`pcaMethods::pca()`) accommodates missing values common in DIA data.

**Differential Abundance:** One-way ANOVA per protein with Tukey HSD post-hoc correction, stratified by sex and age group. Significance threshold: adjusted p < 0.05.

**Pathway Enrichment:** KEGG over-representation analysis (ORA) using `clusterProfiler::enrichKEGG()` with Entrez Gene ID mapping via `AnnotationDbi`.

---

## 📖 Citation

If you use this pipeline, please cite:

> Pandey RS. Brain DIA Proteomics Analysis Pipeline. GitHub: https://github.com/pandeyravi15/dia-proteomics-pipeline

And the underlying tools:
- [clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141)
- [EnhancedVolcano](https://doi.org/10.18129/B9.bioc.EnhancedVolcano)
- [pcaMethods](https://doi.org/10.1093/bioinformatics/btm069)

---

## 🤝 Contributing

Issues and pull requests are welcome. Please open an issue first to discuss major changes.

---

## 📄 License

MIT — see [LICENSE](LICENSE) for details.

---

## 📬 Contact

**Ravi S Pandey**  
GitHub: [@pandeyravi15](https://github.com/pandeyravi15)
