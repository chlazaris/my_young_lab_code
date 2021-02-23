#!/bin/Rscript

# This is to create a scatterplot - takes as input the two variables to be plotted

# Load the required packages
if (!require(ggpubr)) install.packages(ggpubr)

# Load the required libraries
library(ggpubr)

# Get the input arguments
args <- commandArgs(trailingOnly=TRUE)

# Read in the data
data <- read.table(sprintf("%s", args[1]), header=TRUE, stringsAsFactors=F)

# Get the input name
input <- args[1]
output <- gsub(".tsv","",input)

pdf(paste0(output, "_pearson_scatter.pdf"), useDingbats=FALSE, family="sans")
ggscatter(data, x = "gene_length", y = "FPKM",
          add = "reg.line",                               # Add regression line
          conf.int = TRUE,                                # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+ yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE) +
  stat_cor(method = "pearson", label.x = 3, label.y = 4)  # Add correlation coefficient
dev.off()

pdf(paste0(output, "_spearman_scatter.pdf"), useDingbats=FALSE, family="sans")
ggscatter(data, x = "gene_length", y = "FPKM",
          add = "reg.line",                               # Add regression line
          conf.int = TRUE,                                # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+ yscale("log10", .format = TRUE) + xscale("log10", .format = TRUE) +
  stat_cor(method = "spearman", label.x = 3, label.y = 4)  # Add correlation coefficient
dev.off()
