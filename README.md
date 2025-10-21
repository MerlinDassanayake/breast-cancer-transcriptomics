# TCGA Breast Cancer (BRCA) RNA-seq Analysis Pipeline

**Author:** Merlin Dassanayake  
**Language:** R  
**Key Tools:** TCGAbiolinks, DESeq2, clusterProfiler, fgsea, survminer, tidyverse  
**Dataset:** [The Cancer Genome Atlas (TCGA) – BRCA cohort](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)

---

## Overview

This repository presents a complete RNA-seq analysis workflow using TCGA-BRCA data.  
The pipeline performs data acquisition, quality control, differential gene expression, functional enrichment, and survival analysis.  
All scripts are fully automated and can be executed sequentially using a provided Bash script.

The project demonstrates reproducible bioinformatics practices for transcriptomic analysis in cancer research.

---

## Analysis Workflow

| Step | Script | Description |
|------|---------|-------------|
| 1 | `scripts/01_download_tcga.R` | Downloads TCGA-BRCA RNA-seq and clinical data using `TCGAbiolinks`. |
| 2 | `scripts/02_qc_eda.R` | Performs gene filtering, PCA visualization, and top-variable gene heatmap. |
| 3 | `scripts/03_deseq2_deg.R` | Runs differential expression analysis using `DESeq2` and generates volcano plots. |
| 4 | `scripts/04_enrichment.R` | Conducts GO term enrichment and Hallmark GSEA using `clusterProfiler` and `fgsea`. |
| 5 | `scripts/05_survival.R` | Links top DEGs with clinical outcomes and visualizes Kaplan–Meier survival curves. |

All steps can be executed automatically via:

```bash
bash run_analysis.sh
