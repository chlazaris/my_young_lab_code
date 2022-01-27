#!/bin/Rscript

# Load the necessary libraries
library("biomaRt")

# Read the arguments
args <- commandArgs(trailingOnly=TRUE)

# Get the gene symbols
data <- read.table(sprintf("%s", args[1]), header=F)
colnames(data) <- c("gene_name")

# Get the corresponding ENSEMBL IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id","entrezgene_id"),  
      values = data$gene_name, 
      mart = ensembl)

complete_ids <- complete.cases(ids)

