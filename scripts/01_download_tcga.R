# Load libraries in
library(TCGAbiolinks)

# Query: STAR counts for TCGA-BRCA
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

# Download raw files uing GDC api
GDCdownload(query)

# SummarizedExperiment object: counts + metadata
se <- GDCprepare(query)

# Keep primary tumor (TP) and solid tissue normal (NT) samples
# coldata <- colData(se)
# table(coldata$short_letter_code)

# Save the prepared object for downstream analysis
saveRDS(se, 'data_raw_tcga_brca.rds')
