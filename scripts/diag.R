library(SummarizedExperiment)
library(DESeq2)

se <- readRDS("data/raw/tcga_brca_counts_se.rds")

# Check basic statistics
print(paste("Original dimensions:", dim(se)[1], "genes,", dim(se)[2], "samples"))
print("Sample types in the data:")
print(table(colData(se)$sample_type))
print(table(colData(se)$short_letter_code))

# Check if counts are actually zero
counts <- assay(se)
print(paste("Total counts in entire dataset:", sum(counts)))
print(paste("Range of counts:", min(counts), "-", max(counts)))
print(paste("Number of zero-count genes:", sum(rowSums(counts) == 0)))

# Check a few specific genes
print("Counts for first 5 genes, first 5 samples:")
print(counts[1:5, 1:5])