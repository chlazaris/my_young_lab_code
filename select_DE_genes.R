#!/bin/Rscript
library(dplyr)

# Read the input file 
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
data <- read.csv(sprintf("%s", input), header=T)

# Select the relevant columns 
data2 <- data %>% select(EnsemblID, Symbols, log2FoldChange, padj, medianTPMresistant) %>% filter(padj < 0.05 & medianTPMresistant > 1) %>% data.frame()

# Export what you will need
df <- data2 %>% select(EnsemblID, log2FoldChange, padj)
colnames(df) <- c('EnsemblID','log2FoldChange','pvalue')
write.table(df, 'DE_input_to_rank.tsv', row.names=F, col.names=T, sep="\t", quote=F)
