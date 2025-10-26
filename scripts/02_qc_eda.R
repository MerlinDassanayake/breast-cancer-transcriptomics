# Quality control & exploratory analysis
library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# Load in RangedSummarizedExperiment Object
se <- readRDS("data/raw/tcga_brca_counts_se.rds")

# Basic filtering: keep protein-coding genes and drop low counts
se <- se[rowData(se)$gene_type == 'protein_coding', ]
counts <- assay(se)
keep_genes <- rowSums(counts >= 10) >= 10
se <- se[keep_genes, ]

# Use sample_type to filter further
# Primary Tumor = Tumor and Solid Tissue Normal = Normal
colData(se)$group <- ifelse(colData(se)$sample_type == 'Primary Tumor', 'Tumor',
                            'Normal')

# Build DESeqDataSet with design (Tumor vs Normal)
dds_qc <- DESeqDataSet(se, design = ~ group)

# Variance Stabilizing Transform
vsd <- vst(dds_qc, blind = TRUE)

# Save dds and vsd into processed data
saveRDS(dds_qc, 'data/processed/dds_qc.rds')
saveRDS(vsd, 'data/processed/vsd_qc.rds')

# PCA Analysis
pca <- prcomp(t(assay(vsd)))
pc_df <- as_tibble(pca$x[,1:2], rownames = 'sample')  # Extract first two PCs
pc_df$group <- colData(se)$group  # Append group column

# Simple PCA plot with two PCs
ggplot(pc_df, aes(PC1, PC2, color = group)) + 
  geom_point(size = 2.5) + 
  theme_minimal() +
  labs(title = "PCA Plot - Tumor vs Normal Samples")
ggsave("results/figures/pca_samples.png", width = 7, height = 5, dpi = 300)

# Close open devices after ggsave
# Still having issues with RPlot.pdf creation during bash script run
while(!is.null(dev.list())) dev.off()

# Find top 50 genes with most variance and create heatmap
topvar <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topvar, ]
mat <- t(scale(t(mat)))  # Normalization to compare genes
gene_names <- rowData(vsd)$gene_name[topvar]
rownames(mat) <- gene_names


# Subsample colData samples to simplify heatmap
set.seed(123)
sampled_columns <- sample(1:ncol(mat), 100)
mat_sampled <- mat[, sampled_columns]

# Create annotation column for heatmap (what group each sample is from)
annotation_col <- data.frame(
  Condition = colData(vsd)$group[sampled_columns]
)
rownames(annotation_col) <- colnames(mat_sampled)

# Create heatmap
png("results/figures/heatmap_topvar.png", width = 900, height = 1000, res = 150)
pheatmap(
  mat_sampled,
  annotation_col = annotation_col,
  annotation_colors = list(
    Condition = c("Tumor" = "red", "Normal" = "blue")),
  main = "Top 50 Variable Genes",
  color = colorRampPalette(c("navy", "white", "red"))(50),
  legend = TRUE,
  annotation_legend = TRUE,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Low expression", "Mean", "High expression"),
  show_rownames = TRUE,
  show_colnames = FALSE,
)
dev.off()

print("Script 02 Complete!")
