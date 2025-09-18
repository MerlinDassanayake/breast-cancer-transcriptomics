library(TCGAbiolinks)

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

GDCdownload(query)

se <- GDCprepare(query)

saveRDS(se, 'data_raw_tcga_brca.rds')

dim(se) # genes x samples
head(assay(se)) # counts
table(colData(se)$short_letter_code) # Tumor vs Normal counts

