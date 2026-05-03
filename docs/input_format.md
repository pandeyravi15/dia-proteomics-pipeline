# Input File Format Guide

This document describes the required format for input files used in the DIA Proteomics Pipeline.

---

## 1. Protein Abundance Matrix

**File:** CSV format, placed in `data/`

### Format

- **Rows:** One per protein (or protein-peptide combination)
- **Columns:** First column is `ProteinID`; remaining columns are one per sample
- **Values:** Normalized protein abundance (log2-scale recommended)
- **Missing values:** `NA` is acceptable; proteins with >50% missingness per group should be removed upstream

### Example

```
ProteinID,S001,S002,S003,S004,S005
APP|peptide1,12.3,11.8,12.1,13.0,12.5
MAPT|peptide2,10.5,10.2,10.8,9.8,10.1
APOE,9.1,9.4,9.0,8.9,9.2
```

### Notes

- If your data has `Protein|Peptide` format (pipe-delimited), set `has_peptide_suffix <- TRUE` in the config
- If protein IDs are already clean (no peptide suffix), set `has_peptide_suffix <- FALSE`
- Duplicate proteins (after removing peptide suffix) are collapsed by averaging
- Sample IDs in column headers must exactly match the `X` column in the metadata file

---

## 2. Sample Metadata File

**File:** CSV format, placed in `data/`

### Required Columns

| Column | Description | Example Values |
|--------|-------------|----------------|
| `X` | Sample ID (must match column names in abundance matrix) | S001, S002 |
| `Genotype` | Group/genotype label | Control, CaseA, CaseB |
| `Sex` | Biological sex | Male, Female |
| `Age` | Numeric timepoint | 4, 12 |

### Example

```
X,Sex,Genotype,Age
S001,Female,Control,4
S002,Male,Control,4
S003,Female,CaseA,4
S004,Male,CaseA,4
S005,Female,CaseB,12
```

### Notes

- Genotype labels must exactly match what you enter in the `control_genotype` and `case_genotypes` variables in the config block
- Age should be numeric (not character)
- Additional metadata columns are allowed and will be ignored

---

## 3. Normalization Expectations

This pipeline assumes data has already been normalized. Recommended upstream steps:

1. **Missing value filtering:** Remove proteins with >50% missingness per group
2. **Batch correction:** Use TAMPOR, ComBat, or similar tools if data spans multiple batches
3. **Scale:** Log2-transformed abundance values work best for downstream linear models

If your data is not yet normalized, add a normalization chunk before the PCA section.

---

## 4. Species-specific Settings

| Organism | `organism_db` | `kegg_organism` |
|----------|---------------|-----------------|
| Mouse | `"org.Mm.eg.db"` | `"mmu"` |
| Human | `"org.Hs.eg.db"` | `"hsa"` |
| Rat | `"org.Rn.eg.db"` | `"rno"` |
| Zebrafish | `"org.Dr.eg.db"` | `"dre"` |

Install the relevant annotation package via:
```r
BiocManager::install("org.Hs.eg.db")  # for human
```
