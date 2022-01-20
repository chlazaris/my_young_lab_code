#!/bin/Rscript

# Load the required libraries
library(dplyr)

# Combine the diff expression results in a dataframe
filenames = list.files(pattern="*diff_expressed.tsv")

# Combine the differential expression results
# in a single dataframe
files = lapply(filenames, read.delim)
df = do.call(rbind, files)

# For the genes that occur more than one times
# (differentially expressed in more than one condition)
# keep the largest absolute log fold change
result <- df %>%
  group_by(Gene) %>%
  summarise(absLog2FoldChange = max(abs(log2FoldChange), na.rm=TRUE))
result <- data.frame(result[rev(order(result$absLog2FoldChange)),])

# Write the result to a file
write.table(result, "abs_log2FC_sorted_diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)
