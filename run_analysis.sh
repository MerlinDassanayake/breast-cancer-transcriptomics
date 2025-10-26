#!/bin/bash
set -e

echo "Starting TCGA-BRCA pipeline"
#Separate R sessions

echo "Step 1: Downloading TCGA data"
Rscript scripts/01_download_tcga.R

echo "Step 2: Quality Control and exploratory analysis"
Rscript scripts/02_qc_eda.R

echo "Step 3: Differential expression analysis"
Rscript scripts/03_deseq2_deg.R

echo "Step 4: Functional enrichment analysis"
Rscript scripts/04_enrichment.R

echo "Step 5: Survival analysis"
Rscript scripts/05_survival.R

echo "--------------------------------------------"
echo "Pipeline Complete!"
echo "All results saved to 'results/'"
