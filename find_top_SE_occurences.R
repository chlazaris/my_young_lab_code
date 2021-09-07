#!/bin/Rscript

# Load the required libraries
library(dplyr)

data <- read.table("SE_genes_per_patient.tsv", header=F, sep="\t", stringsAsFactors=FALSE)
colnames(data) <- c("donor","se_gene")
data2 <- data %>% count(se_gene, sort=TRUE)

# Write to a file
write.table(data2, "top_SE_gene_occurences.tsv", col.names=F, row.names=F, sep="\t", quote=FALSE)
