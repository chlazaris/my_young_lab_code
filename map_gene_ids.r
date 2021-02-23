#!/bin/Rscript

# This is a script that maps gene IDs 

# Load the required libraries
library(clusterProfiler)
# Get the annotations for human
library(org.Hs.eg.db)
# Get the annotations for mouse
library(org.Mm.eg.db)

# Check for the number of arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!= 1) {print("Please provide correct number of arguments..."); quit(save="no")}

# Read in a data frame with the gene IDs
gene_ID <- read.table(sprintf("%s", args[1]), header=F, sep="\t", stringsAsFactors=F)
colnames(gene_ID) <- c("id")

# Get the IDs
ids <- bitr(gene_ID$id, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(ids) <- c("symbol", "entrezid")

# Write the table with the ids
write.table(ids, "gene_id_conversion_table.tsv", row.names=F, col.names=T, sep="\t", quote=F)
