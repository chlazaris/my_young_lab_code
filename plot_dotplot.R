#!/bin/Rscript
library(ggplot2)
library(forcats) ## for reordering the factor
library(dplyr)

# Get the number of arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("Please provide input .csv file"); q(save="no")}

# Get the input file
input <- args[1]

# Read in the .csv file
data <- read.csv(sprintf("%s", input), header=T, stringsAsFactors=FALSE)
data$GeneRatio <- data$intersection_size/data$term_size 
head(data)

# Plot the data for KEGG
data2 <- subset(data, data$source=="KEGG")

## plot
pdf("KEGG_pathway_enrichment_analysis_with_SE_target_genes.pdf")
ggplot(data2, aes(x = GeneRatio, y = fct_reorder(term_name, GeneRatio))) +
               geom_point(aes(size = GeneRatio, color = adjusted_p_value)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle("Pathway enrichment analysis")
dev.off()
