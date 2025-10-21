# Data acquisition
library(TCGAbiolinks)
library(SummarizedExperiment)

# Query for gene expression
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

# Download raw files using GDC api
GDCdownload(query)

# Save to summarizedexperiment object
se <- GDCprepare(query)

# Keep primary tumor (TP) and solid tissue normal (NT) samples
sample_types <- c("Primary Tumor", "Solid Tissue Normal")
se <- se[, colData(se)$sample_type %in% sample_types]

# Clinical data query
clin_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = 'Clinical',
  data.type = 'Clinical Supplement',
  data.format = 'BCR Biotab'
)

GDCdownload(clin_query)
clin <- GDCprepare(clin_query)

# Save the prepared file objects for downstream analysis
saveRDS(se, file = 'data/raw/tcga_brca_counts_se.rds')
saveRDS(clin, file = 'data/raw/tcga_brca_clinical.rds')

print("Script 01 Complete!")
