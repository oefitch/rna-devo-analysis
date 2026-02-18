# Load necessary libraries

setwd() #ADD YOUR WORKING DIRECTORY

library(sva)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(rtracklayer)

# List of sample IDs
samples <- c() #Add a list of your sample names, for this script the files are titled "SAMPLENAME_RNA_counts.csv"
               #you will need to update line 24 if your files follow a different naming convention
batch <- c() #Add the batch numbers that matches your samples in the same order as you listed samples
tissue <- c() #Add different tissue types, or ignore this if you don't have multiple tissue types
stage <- c() #Add the developmental stages for your samples, sometimes developmental stages are complicated so I gave them a number value
cov_mat <- cbind (tissue, stage) #you can remove this and make group = stage and remove the covar_mod from line 41)

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

# Remove special lines (e.g., "__no_feature", "__ambiguous", etc.) that is a result of HTseq
toRmv <- which(grepl("__", row.names(geneCounts)))
geneCounts <- geneCounts[-toRmv,]

# Assign sample names as column names
colnames(geneCounts) <- samples

# Reduce batch effects 
counts <- geneCounts
adjusted_counts <- ComBat_seq(counts = counts, batch = batch, group = NULL, covar_mod = cov_mat, 
                              full_mod = TRUE, shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL)
