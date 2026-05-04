# Brain Multi-Omics Analysis Pipeline

[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%3E%3D3.16-green)](https://bioconductor.org/)

A reproducible R-based pipeline for end-to-end analysis of **DIA proteomics** and **bulk RNA-seq transcriptomics** data from mouse or human brain tissue. Originally developed for Alzheimer's disease mouse models but generalized for any case-control multi-omics study with multiple genotypes, sexes, and timepoints.

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

### Notebook 02 — Proteomics Human Translation Analysis

| Step | Analysis |
|------|----------|
| Linear Modeling | `abundance ~ Genotype + Age + Sex` per protein with covariate control |
| Module Correlation (logFC) | Pearson correlation of ANOVA fold-changes with human AD co-expression modules |
| Module Correlation (LM) | Pearson correlation of LM estimates with human AD co-expression modules |
| BioDomain ORA | GO over-representation analysis annotated with AMP-AD biological domains |
| BioDomain GSEA | Full ranked GSEA annotated with AMP-AD biological domains |
| Sub-BioDomain | Portrait and landscape plots across AD biological sub-domains |

### Notebook 03 — RNA-seq Transcriptomics Pipeline

| Step | Analysis |
|------|----------|
| Data Inspection | Metadata loading, count matrix QC, transgene merging (optional) |
| Isoform Expression | Visualization of specific gene isoforms e.g. human MAPT (optional) |
| Sex Validation | Xist / Ddx3y expression to confirm sex assignments |
| PCA | DESeq2 VST-based PCA + PCAtools biplot and eigencorplot |
| Differential Expression | DESeq2 pairwise comparisons, stratified by sex and age |
| Visualization | Volcano plots, DEG barplots |
| Functional Enrichment | KEGG pathway ORA via `clusterProfiler` |

### Notebook 04 — RNA-seq Human Translation Analysis

| Step | Analysis |
|------|----------|
| Module Correlation (logFC) | Pearson correlation of DESeq2 logFC with sex-stratified AMP-AD modules |
| Linear Modeling | Per-gene LM on TPM data with genotype, sex, and age covariates |
| Module Correlation (LM) | Pearson correlation of LM estimates with AMP-AD modules |
| AD Subtype Correlation | Correlation with Nikhil et al. and Neff et al. AD subtype signatures |
| BioDomain Correlation (LM) | LM estimates vs AMP-AD BioDomain reference fold-changes |
| BioDomain ORA | GO over-representation annotated with AMP-AD biological domains |
| BioDomain GSEA | Full ranked GSEA annotated with AMP-AD biological domains |

---

## 🗂️ Repository Structure

```
brain-multiomics-pipeline/
│
├── README.md                                        # This file
├── LICENSE                                          # MIT license
├── .gitignore                                       # Excludes data, results, and large files
│
├── notebooks/
│   ├── 01_DIA_Proteomics_Pipeline.Rmd               # Proteomics QC, PCA, ANOVA, KEGG
│   ├── 02_Human_Translation_Analysis.Rmd            # Proteomics LM, modules, BioDomain
│   ├── 03_RNAseq_Pipeline.Rmd                       # RNA-seq QC, PCA, DESeq2, KEGG
│   └── 04_RNAseq_Human_Translation.Rmd              # RNA-seq LM, modules, subtypes, BioDomain
│
├── R/
│   ├── Helper_Functions_Prot.R                      # Proteomics utility functions
│   └── Helper_Functions.R                           # Transcriptomics utility functions
│
├── data/
│   └── example/
│       ├── example_abundance.csv                    # Example protein abundance matrix
│       └── example_traits.csv                       # Example sample metadata
│
├── figures/                                         # Example output figures
│
├── results/                                         # gitignored — outputs saved here locally
│
└── docs/
    └── input_format.md                              # Detailed input file format guide
```

> **Run order:** Notebooks must be run in sequence within each omics type.
> Notebook 02 requires outputs from 01. Notebook 04 requires outputs from 03.
> The proteomics and transcriptomics pipelines are independent of each other.

---

## ⚙️ Requirements

### R Version
R ≥ 4.2.0

### CRAN Packages
```r
install.packages(c(
  "tidyverse", "gt", "UpSetR", "data.table", "broom",
  "ggpubr", "ggh4x", "corrplot", "cowplot", "ggrepel",
  "matrixStats", "VennDiagram", "readxl", "readr"
))
```

### Bioconductor Packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  # Core
  "AnnotationDbi",     # Gene ID mapping
  "clusterProfiler",   # Pathway enrichment (ORA and GSEA)
  "enrichplot",        # Enrichment visualization
  "ComplexHeatmap",    # Heatmap visualization

  # Proteomics specific
  "pcaMethods",        # Probabilistic PCA (handles missing values)

  # Transcriptomics specific
  "DESeq2",            # Differential expression
  "PCAtools",          # Advanced PCA and eigencorplot
  "DEGreport",         # DEG reporting utilities
  "ReactomePA",        # Reactome pathway analysis
  "limma",             # Linear models for omics

  # Organism annotation — install only what you need
  "org.Mm.eg.db",      # Mouse
  "org.Hs.eg.db",      # Human
  "org.Rn.eg.db"       # Rat
))
```

> **Species note:** Only one organism package is needed per study. Set `organism_db` in the USER CONFIGURATION block of each notebook.

### Optional — Human Translation (Notebooks 02 and 04)

External reference datasets required for human comparison analyses:

| Reference | Source | Used in | Used for |
|-----------|--------|---------|----------|
| AMP-AD co-expression modules | AMP-AD Knowledge Portal | 02, 04 | Module logFC and LM correlation |
| AMP-AD sex-stratified modules | AMP-AD Knowledge Portal | 04 | Sex-specific module correlation |
| BioDomain annotations | [Agora / AMP-AD](https://agora.adknowledgeportal.org) | 02, 04 | GO term → BioDomain mapping |
| Domain color table | AMP-AD | 02, 04 | BioDomain plot styling |
| Nikhil et al. AD subtypes | [Rao et al. 2021](https://pubmed.ncbi.nlm.nih.gov/32492070/) | 04 | AD subtype correlation |
| Neff et al. AD subtypes | [Neff et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33523961/) | 04 | AD subtype correlation |

---

## 🚀 Quick Start

### Proteomics Pipeline (Notebooks 01 → 02)

**Step 1:** Prepare your input files (see `docs/input_format.md`):
- Normalized protein abundance matrix (CSV)
- Sample metadata (CSV)

**Step 2:** Configure and run notebook 01:

```r
# Key settings in the USER CONFIGURATION block:
data_dir         <- "../data/"
abundance_file   <- "your_abundance_file.csv"
traits_file      <- "your_metadata_file.csv"
result_dir       <- "../results/"

control_genotype <- "Your_Control_Genotype"
case_genotypes   <- c("Case_1", "Case_2", "Case_3")
age_groups       <- c(4, 12)
organism_db      <- "org.Mm.eg.db"    # swap for human
kegg_organism    <- "mmu"             # "hsa" for human

rmarkdown::render("notebooks/01_DIA_Proteomics_Pipeline.Rmd")
```

**Step 3:** Configure and run notebook 02 (optional):

```r
# Enable only the analyses you have reference data for
run_module_corr    <- TRUE
run_biodomain_ora  <- TRUE
run_biodomain_gsea <- TRUE
run_lm_analysis    <- TRUE

rmarkdown::render("notebooks/02_Human_Translation_Analysis.Rmd")
```

---

### Transcriptomics Pipeline (Notebooks 03 → 04)

**Step 1:** Prepare your input files:
- RSEM gene count matrix (TSV)
- RSEM gene TPM matrix (TSV)
- Sample metadata (CSV)
- Optionally: isoform TPM matrix (TSV)

**Step 2:** Configure and run notebook 03:

```r
# Key settings in the USER CONFIGURATION block:
data_dir           <- "../data/"
gene_counts_file   <- "rsem.merged.gene_counts.tsv"
gene_tpm_file      <- "rsem.merged.gene_tpm.tsv"
metadata_file      <- "metadata_validated.csv"
result_dir         <- "../results/"

control_genotype   <- "Your_Control_Genotype"
case_genotypes     <- c("Case_1", "Case_2", "Case_3")
age_groups         <- c(4, 12)
organism_db        <- "org.Mm.eg.db"
kegg_organism      <- "mmu"

# Transgene handling — set TRUE only if your model has human transgenes
has_transgenes     <- FALSE

rmarkdown::render("notebooks/03_RNAseq_Pipeline.Rmd")
```

**Step 3:** Configure and run notebook 04 (optional):

```r
# Enable only the analyses you have reference data for
run_module_corr_fc <- TRUE
run_module_corr_lm <- TRUE
run_subtype_corr   <- TRUE
run_lm_analysis    <- TRUE
run_biodomain_ora  <- TRUE
run_biodomain_gsea <- TRUE

rmarkdown::render("notebooks/04_RNAseq_Human_Translation.Rmd")
```

---

## 📁 Input File Formats

### Proteomics — Abundance Matrix

```
ProteinID,Sample1,Sample2,Sample3,...
MAPT|peptide1,12.3,11.8,12.1,...
APP|peptide2,10.5,10.2,10.8,...
```

- Rows = proteins; columns = samples
- Values should be pre-normalized (e.g., TAMPOR-corrected, log2-scale)
- ProteinID can be `Protein` or `Protein|Peptide` format

### All Notebooks — Sample Metadata

```
X,Sex,Genotype,Age
Sample1,Female,Control,4
Sample2,Male,Control,4
Sample3,Female,CaseA,4
```

- `X`: Sample ID — must exactly match column names in abundance/count matrix
- `Sex`: Male / Female
- `Genotype`: group labels matching your config variables
- `Age`: numeric timepoint in months

### Transcriptomics — Count Matrix

Standard RSEM TSV output with `gene_id` as first column and one column per sample. See `docs/input_format.md` for full format details and edge cases.

---

## 🧪 Test with Example Data

```r
# For proteomics notebook 01, set in USER CONFIGURATION:
data_dir       <- "../data/example/"
abundance_file <- "example_abundance.csv"
traits_file    <- "example_traits.csv"
```

Example data contains 20 proteins and 12 samples across 2 genotypes, 2 sexes, and 2 ages. Designed to verify end-to-end pipeline execution — enrichment results are not biologically meaningful at this scale.

---

## 📊 Output Files

| Output | Location | Generated by | Description |
|--------|----------|--------------|-------------|
| HTML Reports | `notebooks/` | All notebooks | Full rendered analysis with all plots |
| Processed Proteomics Data | `data/Processed_DIA_Proteomics.RData` | Notebook 01 | Cleaned abundance matrix + traits |
| ANOVA Results | `results/ANOVA_Results_DIA_Proteomics.Rdata` | Notebook 01 | All DEP results by comparison |
| Proteomics Human Translation | `results/Human_Translation_Results.Rdata` | Notebook 02 | LM, module, BioDomain results |
| Processed RNA-seq Data | `data/ProcessedData_RNAseq.RData` | Notebook 03 | Filtered count + TPM matrices |
| DESeq2 Results | `results/DESeq2_Results_RNAseq.RData` | Notebook 03 | All DEG results by comparison |
| LM Results (RNA-seq) | `results/LM_Results_RNAseq.RData` | Notebook 04 | LM estimates by gene and variant |
| ORA BioDomain (RNA-seq) | `results/ORA_BioDomain_RNAseq.RData` | Notebook 04 | GO-ORA annotated with BioDomains |
| GSEA BioDomain (RNA-seq) | `results/GSEA_BioDomain_RNAseq.RData` | Notebook 04 | GO-GSEA annotated with BioDomains |
| RNA-seq Human Translation | `results/Human_Translation_RNAseq.RData` | Notebook 04 | All human correlation results |

---

## 🔬 Methods Summary

**Proteomics QC:** Sample distributions assessed by boxplot after normalization. Proteins with pipe-delimited IDs summarized to protein level by averaging duplicate entries. Probabilistic PCA (`pcaMethods::pca()`) accommodates missing values common in DIA data.

**Differential Abundance (ANOVA):** One-way ANOVA per protein with Tukey HSD post-hoc correction, stratified by sex and age group. Only group pairs present in the data are tested — invalid comparisons are automatically skipped.

**Transcriptomics QC:** Optional transgene merging for models with human gene inserts (configurable via lookup table). Sex validated using Xist and Ddx3y expression. Low-count genes pre-filtered by minimum TPM threshold.

**Differential Expression (DESeq2):** Pairwise DESeq2 comparisons per sex and age stratum using raw integer counts. Genes annotated with symbols and Entrez IDs via `AnnotationDbi`.

**Linear Modeling:** Per-gene/protein linear models estimate the independent contribution of each genotype while controlling for covariates (Age, Sex). Model formula and coefficients are fully user-configurable. Stratified by age group where relevant.

**Pathway Enrichment:** KEGG over-representation analysis using `clusterProfiler::compareCluster()`. Gene set enrichment analysis (GSEA) using `clusterProfiler::gseGO()` on full ranked gene lists.

**Human Translation:** Mouse model fold-changes and LM estimates correlated (Pearson) with human AMP-AD co-expression module signatures (sex-stratified and combined). GO enrichment results mapped to the AMP-AD BioDomain framework across 19 AD biological domains and sub-domains. AD subtype correlations against Nikhil et al. (ROSMAP, Mayo, MSBB) and Neff et al. (MSBB) reference datasets.

---

## 📖 Citation

If you use this pipeline, please cite:

> Pandey RS. Brain Multi-Omics Analysis Pipeline. GitHub: https://github.com/pandeyravi15/brain-multiomics-pipeline

And the underlying tools:
- [DESeq2](https://doi.org/10.1186/s13059-014-0550-8)
- [clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141)
- [EnhancedVolcano](https://doi.org/10.18129/B9.bioc.EnhancedVolcano)
- [pcaMethods](https://doi.org/10.1093/bioinformatics/btm069)
- [PCAtools](https://doi.org/10.18129/B9.bioc.PCAtools)
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
GitHub: [@pandeyravi15](https://github.com/pandeyravi15)
