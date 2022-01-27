#!/bin/Rscript

# This is a script that takes as input an one-column file
# with unique gene names (one per row) and generates a file
# with all combinations of each gene with all the others 
# (network edges), to be used as input to igraph

# Load the required libraries
library(tidyr)
library(dplyr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)
# Check for the number of arguments and quit 
if (length(args)!=1) {
	print("USAGE: Rscript create_edge_matrix INPUT-FILE ")
	q(save="no")
}

# Specify the input
input=args[1]

# Read the input file 
data <- read.table(sprintf("%s", input), header=F, sep="\t", stringsAsFactors=F)
colnames(data) <- c("genes")

# Create a new column which will be a copy of the first one
data$genes2 <- data$genes
# Create a tibble to use with expand
t <- tibble(g1=as.vector(data$genes), g2=as.vector(data$genes2))
df <- t %>% expand(g1, g2) %>% data.frame() 
# Write the dataframe to an output file
write.table(df, "top_CRC_edges.tsv", row.names=F, col.names=T, sep="\t", quote=F)



