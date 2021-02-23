#!/bin/Rscript

# Load the required libraries
library(dplyr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Get the input
input <- args[1]

# Read in the input file
input <- read.table(sprintf("%s", input), header=TRUE, sep="\t", stringsAsFactors=FALSE)
#colnames(input)

# Generate the input 
df <- input %>% select(Description,P.value,N.Geneset.Peak.Genes,N.Geneset.Genes) %>%
	mutate(negative_log10_of_adjusted_p_value=-log(P.value)) %>% rename(term_name=Description, p_value=P.value, intersection_size=N.Geneset.Peak.Genes, term_size=N.Geneset.Genes)

# Read and subset based on the categories of interest
categories <- read.table("GO_categories.tsv", header=T, sep="\t")
df_final <- inner_join(df, categories, by="term_name")
# Write the output into a file to use for plotting the results
write.table(df_final,"GO_input.tsv",col.names=T,row.names=F,sep="\t",quote=F) 
