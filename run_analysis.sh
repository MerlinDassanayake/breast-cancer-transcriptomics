#!/bin/bash

#Separate R sessions
Rscript scripts/01_download_tcga.R
Rscript scripts/02_qc_eda.R
Rscript scripts/03_deseq2_deg.R
Rscript scripts/04_enrichment.R
Rscript scripts/05_survival.R
