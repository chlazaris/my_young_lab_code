#!/bin/Rscript

# Load the GSEA function
source("/lab/solexa_young/lazaris/code/gsea_function.r")
# Load required libraries
library(dplyr)

# Now run GSEA
S4table = read.csv("./Data/journal.pone.0145322.s007.csv", header=TRUE, skip =1) %>%
  filter(Gene.Symbol != "")
gene_list = S4table$DESeq2.Log2.Fold.Change
names(gene_list) = S4table$Gene.Symbol
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

# Specify the GO file
GO_file = "../Annots/c5.bp.v6.2.symbols.gmt"
res = GSEA(gene_list, GO_file, pval = 0.05)
dim(res$Results)
# Plot the results
res$Plot

