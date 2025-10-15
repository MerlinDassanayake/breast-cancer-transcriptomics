# scripts/05_survival.R
library(tidyverse)
library(survival)
library(survminer)

# Load required objects
vsd <- readRDS("data/processed/vsd_qc.rds")
clin <- readRDS("data/raw/tcga_brca_clinical.rds")

# Prepare clinical: try to get days and status robustly
clin2 <- clin %>%
  mutate(days_to_death = as.numeric(days_to_death),
         days_to_last_follow_up = as.numeric(days_to_last_follow_up),
         days = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
         status = ifelse(tolower(vital_status) == "dead", 1, 0)) %>%
  select(case_id, days, status)

# Map sample barcodes to case_id (first 12 characters)
samples <- colnames(vsd)
case_id <- substr(samples, 1, 12)
expr_mat <- assay(vsd)

# Example: test a single gene's expression (MKI67)
gene_symbol <- "MKI67"
# find row in rowData (if present), otherwise map via Annotation
row_ann <- rowData(vsd)
# if rowData has symbol column (common with GDCprepare)
if("external_gene_name" %in% colnames(row_ann)){
  ens <- rownames(vsd)[which(row_ann$external_gene_name == gene_symbol)[1]]
} else {
  stop("No external gene symbol in rowData; map Ensembl IDs to SYMBOL first.")
}
expr <- expr_mat[ens, ]
df_expr <- tibble(sample = samples, case_id = case_id, expr = expr) %>%
  inner_join(clin2, by = "case_id") %>%
  filter(!is.na(days) & !is.na(status))

# cutoff by median
df_expr <- df_expr %>% mutate(group = ifelse(expr > median(expr, na.rm = TRUE), "High", "Low"))

fit <- survfit(Surv(days, status) ~ group, data = df_expr)
p <- ggsurvplot(fit, data = df_expr, pval = TRUE, risk.table = TRUE)
ggsave("results/figures/km_MKI67.png", p$plot, width = 7, height = 6, dpi = 300)

print("Script 05 - Survival Analysis Complete!")