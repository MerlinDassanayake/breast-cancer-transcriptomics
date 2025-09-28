# Loading in Libraries
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
sample.types <- c("Primary Tumor", "Solid Tissue Normal")
se <- se[, colData(se)$sample_type %in% sample.types]

# Clinical data query
clin.query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = 'Clinical',
  data.type = 'Clinical Supplement',
  data.format = 'BCR Biotab'
)

GDCdownload(clin.query)
clin <- GDCprepare(clin.query)

# Save the prepared file objects for downstream analysis
saveRDS(se, file = 'data/raw/tcga_brca_counts_se.rds')
saveRDS(clin, file = 'data/raw/tcga_brca_clinical.rds')
