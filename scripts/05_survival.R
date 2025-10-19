# Survival analysis
library(tidyverse)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(patchwork)

# Load required objects
vsd <- readRDS("data/processed/vsd_qc.rds")
clin <- readRDS("data/raw/tcga_brca_clinical.rds")
deg_results <- read.csv("results/tables/DEG_tumor_vs_normal_full.csv")

# Get top 5 upregulated genes with lowest p-values
top_up_genes <- deg_results %>%  # Significant and strongly upregulated
  filter(padj < 0.05 & log2FoldChange > 2) %>%
  arrange(padj) %>%
  head(5) %>%
  pull(SYMBOL)

# Get top 5 downregulated genes with lowest p-values  
top_down_genes <- deg_results %>%  # Significant and strongly downregulated
  filter(padj < 0.05 & log2FoldChange < -2) %>%
  arrange(padj) %>%
  head(5) %>%
  pull(SYMBOL)


candidate_genes <- c(top_up_genes, top_down_genes)
# cat("Candidate genes for survival analysis:\n")
# print(candidate_genes)

# Prepare clinical: try to get days and status robustly
clin_clean <- clin$clinical_patient_brca %>%
  mutate(
    days_to_death = suppressWarnings(as.numeric(as.character(death_days_to))),
    days_to_last_follow_up = suppressWarnings(as.numeric(as.character(last_contact_days_to))),
    # Survival time using death days or follow up if not available
    days = ifelse(!is.na(days_to_death) & days_to_death > 0, days_to_death,
                  days_to_last_follow_up),
    status = ifelse(tolower(vital_status) == 'dead', 1, 0)
  ) %>%
  dplyr::select(case_id = bcr_patient_barcode, days, status) %>%
  filter(!is.na(days) & days > 0 & !is.na(status))

# Map sample barcodes to case_id (first 12 characters)
samples <- colnames(vsd)
case_ids <- substr(samples, 1, 12)

# Create expression matrix
expr_mat <- assay(vsd)
row_ann <- rowData(vsd)

# Function to create survival plots for each gene
create_survival_plot <- function(gene_symbol, expr_mat, row_ann, case_ids, clin_clean,
                                 samples, regulated) {
  if ("gene_name" %in% colnames(row_ann)) {
    ens <- rownames(vsd)[which(row_ann$gene_name == gene_symbol)[1]]
  }
  
  expr <- expr_mat[ens, ]
  df_expr <- tibble(sample = samples, case_id = case_ids, expr = expr) %>%
    inner_join(clin_clean, by = "case_id") %>%
    filter(!is.na(days) & !is.na(status))
  
  df_expr <- df_expr %>%
    mutate(group = ifelse(expr > median(expr, na.rm = TRUE), "High", "Low"))
  
  # Kaplan-Meier
  fit <- survfit(Surv(days, status) ~ group, data = df_expr)
  
  # Color based on regulation
  if (regulated == "up") {
    plot_colors <- c("#E41A1C", "#FBB4AE")
  } else {
    plot_colors <- c("#377EB8", "#B3CDE3")
  }
  
  # Create survival plot
  s_plot <- ggsurvplot(fit, data = df_expr, pval = TRUE, title = gene_symbol,
                       palette = plot_colors, legend = 'none')
  
  return(s_plot$plot)
}

up_plots <- list()
down_plots <- list()

for (gene in top_up_genes) {
  plot <- create_survival_plot(gene, expr_mat, row_ann, case_ids, clin_clean, samples,
                               "up")
  if (!is.null(plot)) {
    up_plots[[gene]] <- plot
  }
}

for (gene in top_down_genes) {
  plot <- create_survival_plot(gene, expr_mat, row_ann, case_ids, clin_clean, samples,
                               "down")
  if(!is.null(plot)) {
    down_plots[[gene]] <- plot
  }
}

up_panel <- wrap_plots(up_plots, ncol = 5)
down_panel <- wrap_plots(down_plots, ncol = 5)

final_plot <- up_panel / down_panel +
  plot_annotation(
    title = "Survival Analysis of Top Differentially Expressed Genes in Breast Cancer",
    subtitle = "TCGA BRCA Dataset - Separated by Regulation Direction",
    theme_minimal()
  )

ggsave("results/figures/survival_faceted_final.png", 
       final_plot, 
       width = 16, 
       height = 10, 
       dpi = 300,
       bg = "white")

print("Script 05 - Survival Analysis Complete!")