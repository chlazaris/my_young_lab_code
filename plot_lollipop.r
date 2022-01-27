#!/bin/Rscript

# Load the required libraries
library(ggplot2)
library(forcats) ## for reordering the factor
library(dplyr)

# Get the number of arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("Please provide input .csv file"); q(save="no")}
input <- args[1]

# Get the input
data <- read.csv(sprintf("%s", input), header=T, stringsAsFactors=F)
#data$intersection_size <- as.numeric(as.character(data$intersection_size))
#data$term_size <- as.numeric(as.character(data$term_size))
data$GeneRatio <- as.numeric(as.character(data$intersection_size/data$term_size)) 
head(data)

# Name the output
output1 <- gsub(".csv", "_GeneRatio_sorted_plot.pdf", input)

# Plot based on GeneRatio
pdf(sprintf("%s", output1), width=10, height=5, useDingbats=FALSE)
ggplot(data, aes(x = GeneRatio, y = fct_reorder(term_name, GeneRatio))) +
               geom_point(aes(size = GeneRatio, color = negative_log10_of_adjusted_p_value)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(high="red", low="blue") +
        ylab("GO term") +
        ggtitle("Pathway enrichment")
dev.off()

# Name the output
output2 <- gsub(".csv", "_pval_sorted_plot.pdf", input)

# Plot based on p-value 
pdf(sprintf("%s", output2), width=10, height=5, useDingbats=FALSE)
ggplot(data, aes(x = GeneRatio, y = fct_reorder(term_name, negative_log10_of_adjusted_p_value))) +
               geom_point(aes(size = GeneRatio, color = negative_log10_of_adjusted_p_value)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(high="red", low="blue") +
        ylab("GO term") +
        ggtitle("Pathway enrichment")
dev.off()
