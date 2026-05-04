# Brain DIA Proteomics Analysis Pipeline

[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%3E%3D3.16-green)](https://bioconductor.org/)

A reproducible R-based pipeline for end-to-end analysis of **Data-Independent Acquisition (DIA) proteomics data** from mouse or human brain tissue. Originally developed for Alzheimer's disease mouse models but generalized for any case-control DIA proteomics study with multiple genotypes, sexes, and timepoints.

---

## 📋 What This Pipeline Does

### Notebook 01 — DIA Proteomics Pipeline

| Step | Analysis |
|------|----------|
| QC | Protein abundance matrix loading, boxplot QC, summarization to protein level |
| PCA | Probabilistic PCA (handles missing values) with metadata overlay |
| Differential Abundance | ANOVA + Tukey HSD post-hoc, stratified by sex and age |
| Visualization | Volcano plots, DEP barplots, UpSet overlap plots, expression boxplots |
| Functional Enrichment | KEGG pathway ORA via `clusterProfiler` |

### Notebook 02 — Human Translation Analysis

| Step | Analysis |
|------|----------|
| Linear Modeling | `abundance ~ Genotype + Age + Sex` per protein with covariate control |
| Module Correlation | Pearson correlation of mouse fold-changes with human AD co-expression modules |
| LM × Module | Pearson correlation of LM estimates with human AD co-expression modules |
| BioDomain ORA | GO over-representation analysis annotated with AMP-AD biological domains |
| BioDomain GSEA | Full ranked GSEA annotated with AMP-AD biological domains |
| Sub-BioDomain | Protein-level portrait plots across AD biological sub-domains |

---

## 🗂️ Repository Structure

```
dia-proteomics-pipeline/
│
├── README.md                                   # This file
├── LICENSE                                     # MIT license
├── .gitignore                                  # Excludes data, results, and large files
│
├── notebooks/
│   ├── 01_DIA_Proteomics_Pipeline.Rmd          # QC, PCA, ANOVA, KEGG
│   └── 02_Human_Translation_Analysis.Rmd       # LM, module correlation, BioDomain
│
├── R/
│   └── Helper_Functions_Prot.R                 # All shared utility functions
│
├── data/
│   └── example/
│       ├── example_abundance.csv               # Example protein abundance matrix
│       └── example_traits.csv                  # Example sample metadata
│
├── figures/                                    # Example output figures
│   ├── example_PCA.png
│   ├── example_volcano.png
│   └── example_KEGG.png
│
├── results/                                    # gitignored — output saved here locally
│
└── docs/
    └── input_format.md                         # Detailed input file format guide
```

---

## ⚙️ Requirements

### R Version
R ≥ 4.2.0

### CRAN Packages
```r
install.packages(c(
  "tidyverse", "gt", "UpSetR", "EnhancedVolcano", "data.table",
  "ggpubr", "ggh4x", "corrplot", "cowplot", "ggrepel"
))
```

### Bioconductor Packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "pcaMethods",        # Probabilistic PCA
  "clusterProfiler",   # Pathway enrichment (ORA and GSEA)
  "AnnotationDbi",     # Gene ID mapping
  "enrichplot",        # Enrichment visualization
  "ComplexHeatmap",    # Heatmap visualization
  "org.Mm.eg.db",      # Mouse gene annotations
  "org.Hs.eg.db"       # Human gene annotations (swap for human studies)
))
```

> **Species note:** Install `org.Mm.eg.db` for mouse, `org.Hs.eg.db` for human, or `org.Rn.eg.db` for rat. Only one is needed per study — set it in the USER CONFIGURATION block.

### Optional — Human Translation Analysis (Notebook 02)

The following external reference datasets are required only if running notebook 02:

| Reference | Source | Used for |
|-----------|--------|----------|
| Co-expression module signatures | Emory TMT proteomics / AMP-AD | Module correlation |
| BioDomain annotations | [Agora / AMP-AD](https://agora.adknowledgeportal.org) | BioDomain ORA and GSEA |
| Domain color table | AMP-AD | BioDomain plot styling |

---

## 🚀 Quick Start

### Step 1: Prepare Your Input Files

You need two CSV files (see `docs/input_format.md` for details):

**1. Protein Abundance Matrix**
- Rows = proteins, Columns = samples
- Values should be pre-normalized (e.g., median-normalized, TAMPOR-corrected)
- ProteinID column can be `Protein` or `Protein|Peptide` format

**2. Sample Metadata**
- One row per sample
- Required columns: Sample ID, Genotype, Sex, Age

### Step 2: Configure Notebook 01

Open `notebooks/01_DIA_Proteomics_Pipeline.Rmd` and edit the **USER CONFIGURATION** block:

```r
# ── PATHS ──────────────────────────────────────────
data_dir         <- "../data/"
abundance_file   <- "your_abundance_file.csv"
traits_file      <- "your_metadata_file.csv"
result_dir       <- "../results/"

# ── STUDY DESIGN ───────────────────────────────────
control_genotype <- "Your_Control_Genotype"
case_genotypes   <- c("Case_1", "Case_2", "Case_3")
age_groups       <- c(4, 12)

# ── SPECIES ────────────────────────────────────────
organism_db      <- "org.Mm.eg.db"   # mouse
kegg_organism    <- "mmu"            # "hsa" for human
```

### Step 3: Run Notebook 01

In RStudio, click **Knit** or run:

```r
rmarkdown::render("notebooks/01_DIA_Proteomics_Pipeline.Rmd")
```

This generates:
- `results/ANOVA_Results_DIA_Proteomics.Rdata`
- `data/Processed_DIA_Proteomics.RData`
- An HTML report with all QC, DEA, and KEGG results

### Step 4: Configure and Run Notebook 02 (Optional)

Open `notebooks/02_Human_Translation_Analysis.Rmd` and set the analysis switches:

```r
# Enable only the analyses you have reference data for
run_module_corr    <- TRUE    # requires module_data_file
run_biodomain_ora  <- TRUE    # requires biodom_path + domain_colors_path
run_biodomain_gsea <- TRUE    # requires biodom_path + domain_colors_path
run_lm_analysis    <- TRUE    # requires no external data

# Provide paths to human reference data
module_data_file   <- "path/to/module_reference.RData"
biodom_path        <- "path/to/annotated_biodomains.rds"
domain_colors_path <- "path/to/domain_colors.csv"
```

Then run:

```r
rmarkdown::render("notebooks/02_Human_Translation_Analysis.Rmd")
```

> **Note:** Notebook 02 requires output from notebook 01. Always run notebook 01 first.

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

- `X`: Sample ID — must exactly match column names in the abundance matrix
- `Sex`: Male / Female
- `Genotype`: group labels for each sample
- `Age`: numeric timepoint

See `docs/input_format.md` for full format details, edge cases, and species-specific settings.

---

## 🧪 Test with Example Data

```r
# In the notebook 01 USER CONFIGURATION block, set:
data_dir       <- "../data/example/"
abundance_file <- "example_abundance.csv"
traits_file    <- "example_traits.csv"
```

Example data contains 20 proteins and 12 samples across 2 genotypes, 2 sexes, and 2 ages. It is designed to verify the pipeline runs end-to-end — enrichment results will not be biologically meaningful at this scale.

---

## 📊 Output Files

| Output | Location | Generated by | Description |
|--------|----------|--------------|-------------|
| HTML Report (QC + DEA) | `notebooks/` | Notebook 01 | Full rendered analysis with all plots |
| Processed Data | `data/Processed_DIA_Proteomics.RData` | Notebook 01 | Cleaned abundance matrix + traits |
| ANOVA Results | `results/ANOVA_Results_DIA_Proteomics.Rdata` | Notebook 01 | All DEP results by model |
| HTML Report (Human) | `notebooks/` | Notebook 02 | Module correlation, BioDomain, LM results |
| Human Translation Results | `results/Human_Translation_Results.Rdata` | Notebook 02 | LM estimates, module/BioDomain correlations |

---

## 🔬 Methods Summary

**Quality Control:** Sample distributions are assessed by boxplot after normalization. Proteins with pipe-delimited IDs are summarized to protein level by averaging duplicate entries. Probabilistic PCA (`pcaMethods::pca()`) handles missing values common in DIA data.

**Differential Abundance (ANOVA):** One-way ANOVA per protein with Tukey HSD post-hoc correction, stratified by sex and age group. Only group pairs present in the data are tested — invalid comparisons are automatically skipped.

**Linear Modeling:** Per-protein linear models (`abundance ~ Genotype + Age + Sex`) estimate the independent contribution of each genotype while controlling for covariates. Model formula and coefficients are fully user-configurable.

**Pathway Enrichment:** KEGG over-representation analysis (ORA) using `clusterProfiler::enrichKEGG()`. Entrez Gene IDs are mapped via `AnnotationDbi`.

**Human Translation:** Mouse model fold-changes and LM estimates are correlated (Pearson) with human AD brain proteomics co-expression module signatures. GO enrichment results are mapped to the AMP-AD BioDomain framework to identify disease-relevant biological processes.

**BioDomain GSEA:** Gene set enrichment analysis using `clusterProfiler::gseGO()` on the full ranked protein list, with GO terms annotated to AD biological domains and sub-domains.

---

## 📖 Citation

If you use this pipeline, please cite:

> Pandey RS. Brain DIA Proteomics Analysis Pipeline. GitHub: https://github.com/[your-username]/dia-proteomics-pipeline

And the underlying tools:
- [clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141)
- [EnhancedVolcano](https://doi.org/10.18129/B9.bioc.EnhancedVolcano)
- [pcaMethods](https://doi.org/10.1093/bioinformatics/btm069)
- [AMP-AD BioDomain Framework](https://agora.adknowledgeportal.org)

---

## 🤝 Contributing

Issues and pull requests are welcome. Please open an issue first to discuss major changes.

---

## 📄 License

MIT — see [LICENSE](LICENSE) for details.

---

## 📬 Contact

**Ravi S Pandey**  
GitHub: [@your-username](https://github.com/your-username)
