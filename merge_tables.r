#!/bin/Rscript

# Load the required libraries
library(dplyr)

# Get the unique genes regulated by TEs
te_genes <- read.table("TE_regulated_genes_uniq.tsv", header=F, sep="\t", stringsAsFactors=F)
colnames(te_genes) <- c("refseq", "geneSymbol")
te_genes$category <- c("TE")

se_genes <- read.table("SE_regulated_genes_uniq.tsv", header=F, sep="\t", stringsAsFactors=F)
colnames(se_genes) <- c("refseq", "geneSymbol")
se_genes$category <- c("SE")

# Now we have all the genes that are regulated either by SE or TE based on Sabari et al.
all_genes <- rbind(se_genes, te_genes)
head(all_genes)

# Read the gene mapping
mapping <- read.table("mm9_gene_mapping_complete.tsv", header=T, sep="\t", stringsAsFactors=F)

# Combine with the gene mapping that has NM identifiers to get the Geneid
genes_annot <- inner_join(all_genes, mapping, by='refseq')
head(genes_annot)
dim(genes_annot)

