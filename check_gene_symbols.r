#!/bin/Rscript

# This is a script that checks gene symbols against the latest version of
# HGNC database and gives suggestions for how they could be corrected

# Load the required libraries
if (!require(HGNChelper)) install.packages(HGNChelper)
library(dplyr)

# Read in the data
data <- read.table("7878_7879_common_target_genes.tsv", header=F, sep="\t")
colnames(data) <- c("genes")

# Extract the symbols
symbols <- data$genes

# Get the latest version of gene symbols
new_human_map <- getCurrentHumanMap()

# Check the gene symbols based on the latest map
df <- checkGeneSymbols(symbols, map=new_human_map, species="human")
# Get only the ones that have suggested symbols
df2 <- df[!is.na(df$Suggested.Symbol),]

# Extract the columns that you want
df3 <- df2 %>% select(x, Suggested.Symbol) %>% data.frame()
colnames(df3) <- c("geneSymbol", "suggested_geneSymbol")
write.table(df3, "geneSymbol_mapping.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)
