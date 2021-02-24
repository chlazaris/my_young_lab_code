#!/bin/Rscript

# Load the required libraries
library(stringr)
library(ggplot2)

# Function to test whether a certain gene
# is a member of a CRC or not
search_for_factor <- function(x){
	x_with_no_brackets <- gsub("\\[|\\]","",x)
	# Convert it
	x_vector <- unlist(str_split(x_with_no_brackets, ', '))
	result <- ifelse(gene %in% x_vector==TRUE, 1, 0)
}

# Get the command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Get the input
input_file <- args[1]
gene <- as.vector(as.character(args[2]))

# Read input
data <- read.table(sprintf("%s", args[1]), header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Name the columns
colnames(data) <- c("genes","CRC_score","network_size")

# Convert the column with the members of each
# CRC to a dataframe
network_genes <- data.frame(data$genes)
colnames(network_genes) <- c("net_genes")

# Check if the gene is present in the corresponding network
present_in_network <- apply(network_genes, 1, search_for_factor)

# Add the new column with network presence to the data
data$gene_present <- as.factor(present_in_network)
# Add a column with the CRC network rank based on CRC score
data$network_rank <- c(1:nrow(data))
data$network_rank <- sort(data$network_rank, decreasing=FALSE)

# Specify the name of the outfile
outfile <- paste0(gene, "_CRCmapper_results.pdf")

# Plot the data
pdf(outfile, useDingbats=FALSE)
ggplot(data, aes(x=network_rank, y=CRC_score)) +
    geom_segment(aes(xend=network_rank, yend=0), alpha=0.1) +
    geom_point(aes(color=gene_present)) +
    scale_color_manual(labels=c("absent", "present"), values=c("#d3d3d3", "#FF0000")) +
    theme_bw() +
    xlab("Network rank") +
    ylab("CRCmapper score")+ 
    guides(color=guide_legend(title=gene, reverse=TRUE)) +
    scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) 
dev.off()
