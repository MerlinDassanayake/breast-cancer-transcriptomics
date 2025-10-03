library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)

dds <- readRDS('data/processed/dds.rds')
vsd <- readRDS('data/processed/vsd.rds')

# If you have covariates like batch/age, update design before running DESeq(dds)
# e.g., design(dds) <- ~ batch + group
design(dds) <- ~ group

dds <- DESeq(dds)

# Use apeglm (or ashr) for shrinkage (more stable log2FC)
res_shr <- lfcShrink(dds, contrast = c('group', 'Tumor', 'Normal'),
                     type = 'apeglm')

# Make tidy results table
res_df <- as.data.frame(res_shr) %>%
  rownames_to_column(var = 'ensembl') %>%
  mutate(ENSEMBL = gsub("\\..*", "", ensembl)) # Remove Ensembl version suffix if present

# Map to gene symbols and Entrez IDs
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = unique(res_df$ENSEMBL),
                             columns = c("SYMBOL","ENTREZID"),
                             keytype = "ENSEMBL")

res_df <- left_join(res_df, map, by = c("ENSEMBL" = "ENSEMBL"))

# Save full table
write.csv(res_df %>% arrange(padj), "results/tables/DEG_tumor_vs_normal_full.csv", row.names = FALSE)

# Quick volcano plot (label top hits)
res_plot <- res_df %>% filter(!is.na(padj))
top_hits <- res_plot %>% arrange(padj) %>% slice_head(n=10)

p <- ggplot(res_plot, aes(log2FoldChange, -log10(padj))) + 
  geom_point(alpha=0.4) + 
  geom_vline(xintercept = c("-1,1"), linetype="dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  ggrepel::geom_text_repel(data = top_hits, aes(label = SYMBOL), size = 2) +
  theme_minimal()
ggsave("results/figures/volcano.png", p, width = 7, height = 5, dpi = 300)