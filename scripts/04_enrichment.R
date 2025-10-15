# Functional enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(tidyverse)

res_df <- read.csv("results/tables/DEG_tumor_vs_normal_full.csv")

# ORA: use ENTREZ IDs of up/down significant genes
sig_up <- res_df %>% filter(!is.na(ENTREZID) & padj < 0.05 & log2FoldChange > 1) %>%
  pull(ENTREZID) %>% 
  unique()


sig_down <- res_df %>% filter(!is.na(ENTREZID) & padj < 0.05 & log2FoldChange < -1) %>%
  pull(ENTREZID) %>% 
  unique()

ego_up <- enrichGO(gene = sig_up, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
ego_dn <- enrichGO(gene = sig_down, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)

write.csv(as.data.frame(ego_up), "results/tables/GO_up.csv", row.names = FALSE)
write.csv(as.data.frame(ego_dn), "results/tables/GO_dn.csv", row.names = FALSE)

# GSEA using MSigDB Hallmarks (rank all genes by a statistic)
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Build a ranking: sign(logFC) * -log10(padj) (avoid padj==0 by adding tiny eps)
res_rnk <- res_df %>% filter(!is.na(SYMBOL)) %>%
  mutate(stat = sign(log2FoldChange) * -log10(pmax(padj, 1e-300))) %>%
  group_by(SYMBOL) %>% summarize(stat = max(stat, na.rm = TRUE)) %>%
  arrange(desc(stat))

ranks <- deframe(res_rnk) # named numeric vector: names = SYMBOL, values = stat
fgsea_res <- fgsea(pathways = msig_hallmark, stats = ranks, minSize = 15, maxSize = 500)
fgsea_res_clean <- fgsea_res %>%
  mutate(leadingEdge = map_chr(leadingEdge, ~paste(., collapse = ',')))
fgsea_res_clean <- fgsea_res_clean %>% arrange(padj)
write.csv(fgsea_res_clean, "results/tables/gsea_hallmark.csv", row.names = FALSE)

# Plot top pathway using clean dataframe
top_path <- fgsea_res_clean$pathway[1]
png("results/figures/gsea_top.png", width = 800, height = 600, res = 150)
p <- plotEnrichment(msig_hallmark[[top_path]], ranks) + 
  ggtitle(top_path) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
print(p)
dev.off()

print("Script 04 - Functional Enrichment Complete!")