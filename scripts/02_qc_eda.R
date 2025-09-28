library(SummarizedExperiment)
library(DESeq2); library(tidyverse)
library(ComplexHeatmap)

se <- readRDS("data/raw/tcga_brca_counts_se.rds")

# Basic filtering: keep protein-coding genes and drop low counts
if("gene_type" %in% colnames(rowData(se))){
  se <- se[rowData(se)$gene_type == "protein_coding", ]
}
counts <- assay(se)
keep_genes <- rowSums(counts >= 10) >= 10
se <- se[keep_genes, ]

# Use sample_type to filter further
colData(se)$group <- ifelse(coldata$sample_type == 'Primary Tumor', 'Tumor', 'Normal')

# Build DESeqDataSet with design (Tumor vs Normal)
coldata <- as.data.frame(colData(se))

dds <- DESeqDataSet(se, design = ~ group)

# Size factors + variance stabilizing transform
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

# Save dds and vsd into processed data
saveRDS(dds, 'data/processed/dds.rds')
saveRDS(vsd, 'data/processed/vsd.rds')

# Sample QC: PCA
pca <- prcomp(t(assay(vsd)))
pc_df <- as.tibble(pca$x[,1:2], rownames = 'sample')
pc_df$group <- coldata$group

ggplot(pc_df, aes(PC1, PC2, color = group)) + 
  geom_point(size = 2.5) + 
  theme_minimal() +
  labs(title = "PCA Plot - Tumor vs Normal Samples")
ggsave("results/figures/pca_samples.png", width = 7, height = 5, dpi = 300)

# Top-variable genes heatmap (quick look)
topvar <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 500)
mat <- assay(vsd)[topvar, ]
mat <- t(scale(t(mat)))

# FIX: Corrected variable name from map to mat
png("results/figures/heatmap_topvar.png", width = 1000, height = 1000, res = 150)
draw(Heatmap(mat, show_row_names = FALSE, name = 'zscore'))
dev.off()

print("Analysis completed successfully!")