# Quality control & exploratory analysis
library(SummarizedExperiment)
library(DESeq2); library(tidyverse)
library(ggplot2)
library(pheatmap)

# Load in RangedSummarizedExperiment Object
se <- readRDS("data/raw/tcga_brca_counts_se.rds")

# Basic filtering: keep protein-coding genes and drop low counts
if("gene_type" %in% colnames(rowData(se))){
  se <- se[rowData(se)$gene_type == "protein_coding", ]
}
counts <- assay(se)
keep_genes <- rowSums(counts >= 10) >= 10
se <- se[keep_genes, ]

# Use sample_type to filter further
# Primary Tumor = Tumor and Solid Tissue Normal = Normal
colData(se)$group <- ifelse(colData(se)$sample_type == 'Primary Tumor', 'Tumor', 'Normal')

# Build DESeqDataSet with design (Tumor vs Normal)
# coldata <- as.data.frame(colData(se))

dds_qc <- DESeqDataSet(se, design = ~ group)

vsd <- vst(dds_qc, blind = TRUE)

# Save dds and vsd into processed data
saveRDS(dds_qc, 'data/processed/dds_qc.rds')
saveRDS(vsd, 'data/processed/vsd_qc.rds')

# PCA Analysis
pca <- prcomp(t(assay(vsd)))
pc_df <- as_tibble(pca$x[,1:2], rownames = 'sample')
pc_df$group <- colData(se)$group

# Simple PCA plot with two PC's
ggplot(pc_df, aes(PC1, PC2, color = group)) + 
  geom_point(size = 2.5) + 
  theme_minimal() +
  labs(title = "PCA Plot - Tumor vs Normal Samples")
ggsave("results/figures/pca_samples.png", width = 7, height = 5, dpi = 300)

# Find top 50 variable genes and create heatmap
topvar <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topvar, ]
mat <- t(scale(t(mat)))
gene_names <- rowData(vsd)$gene_name[topvar]
rownames(mat) <- gene_names

# Subsample colData samples to simplify heatmap
set.seed(123)
sampled_columns <- sample(1:ncol(mat), 100)
mat_sampled <- mat[, sampled_columns]

png("results/figures/heatmap_topvar.png", width = 900, height = 1000, res = 150)
pheatmap(mat_sampled,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = as.data.frame(colData(vsd)[sampled_columns, "group", drop = FALSE]))
dev.off()

print("Script 02 Complete!")
