# Load necessary libraries

setwd("/mnt/research/FishEvoDevoGeno/Fitch/tail_multiomics/RNA_analysis/batch_effects")

library(sva)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(rtracklayer)

# List of sample IDs
samples <- c("CF1_1", "CF2", "CF3", "CF4", "CF5", "CF6", "CF7", "CF8", "CF9", 
             "CR1", "CR2", "CR3", "CR4", "CR5", "CR6", "CR7", "CR8", "CR9", 
             "NE1_1", "NE2", "NE3", "NE4", "NE5", "NE6", "NE7", "NE8", "NE9")

# Initialize an empty list to store data frames
counts_list <- list()

# Read all count files into a list
for (sample in samples) {
  filename <- paste0(sample, "_RNA_counts.csv")
  counts_list[[sample]] <- read.table(filename, sep = "\t", as.is = TRUE, header = FALSE)
}

# Combine into a single data frame
geneCounts <- do.call(cbind, lapply(counts_list, function(df) df[,2]))  # Extracts count column
row.names(geneCounts) <- counts_list[[1]][,1]  # Use gene names from the first dataset

# Remove special lines (e.g., "__no_feature", "__ambiguous", etc.)
toRmv <- which(grepl("__", row.names(geneCounts)))
geneCounts <- geneCounts[-toRmv,]

# Assign sample names as column names
colnames(geneCounts) <- samples

# Reduce batch effects 
counts <- geneCounts
batch <- c(2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2)
adjusted_counts <- ComBat_seq(counts = counts, batch = batch, group = NULL, covar_mod = NULL, full_mod = TRUE, shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL)

# Create a condition vector based on prefixes
condition <- ifelse(grepl("^CF", samples), "CF",
             ifelse(grepl("^CR", samples), "CR", "NE"))

# Convert condition to a factor
condition <- factor(condition, levels = c("CF", "CR", "NE"))

# Create DGEList object

dge <- DGEList(counts = adjusted_counts, group = condition)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize data using TMM normalization
dge <- calcNormFactors(dge)

# Generate design matrix for differential expression analysis

design <- model.matrix(~condition)

# Estimate dispersions

dge <- estimateDisp(dge, design)

# Perform differential expression analysis using exact test
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit)

# Extract differentially expressed genes
results <- topTags(qlf, n = Inf)$table

# Save results to a CSV file
write.csv(results, file = "edgeR_DE_results.csv")

# Print summary of DE results
summary(decideTests(qlf))

# Plot MDS (Multidimensional Scaling) for sample similarity
plotMDS(dge, col = as.numeric(condition))
