# Differential expression analysis
library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)

# Read in pre-processed qc DESeq2 dataset object
dds_qc <- readRDS("data/processed/dds_qc.rds")

# Run full DESeq2 analysis pipeline
dds <- DESeq(dds_qc)

# Check available coef
print(resultsNames(dds))

# Extract results with apeglm shrinkage
res_shr <- lfcShrink(dds, coef = "group_Tumor_vs_Normal", type = "apeglm")

# Make a tidy results table
res_df <- as.data.frame(res_shr) %>%
  rownames_to_column(var = "ensembl_with_version") %>%
  mutate(ENSEMBL = gsub("\\..*", "", ensembl_with_version))  # remove Ensembl version suffix if present

# Map to gene symbols and Entrez IDs
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = unique(res_df$ENSEMBL),
                             columns = c("SYMBOL","ENTREZID"),
                             keytype = "ENSEMBL")

res_df <- left_join(res_df, map, by = "ENSEMBL") %>% 
  select(-ensembl_with_version)

# Save full table
write.csv(res_df %>% arrange(padj), "results/tables/DEG_tumor_vs_normal_full.csv", row.names = FALSE)

# Quick volcano plot (label top hits)
res_plot <- res_df %>% filter(!is.na(padj))
top_hits <- res_plot %>% arrange(padj) %>% slice_head(n = 10)

p <- ggplot(res_plot, aes(log2FoldChange, -log10(padj))) +
  geom_point(alpha = 0.4) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(data = top_hits, aes(label = SYMBOL), size = 2) +
  theme_minimal()
ggsave("results/figures/volcano.png", p, width = 7, height = 5, dpi = 300)

print("Script 03 - DEG Analysis Complete!")