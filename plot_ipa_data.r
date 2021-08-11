#!/bin/Rscript

# Load the required packages and libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# Get the arguments
input <- commandArgs(trailingOnly=TRUE)
print(input)

# Read in the data
data <- read.table(sprintf("%s", input), sep=",", header=TRUE, check.name=FALSE)

# Check the data
head(data)

# Calculate -log10pvalue and add it to the dataset
data2 <- data %>% mutate(minuslog10pval=-log10(p_value))
head(data2)

pdf("top_upstream_regulators.pdf", useDingbats=FALSE)
ggbarplot(data2, x="Upstream Regulator", y="Activation z-score", fill="minuslog10pval", sort.val = "desc", legend.title = "-log10 p-value", sort.by.groups = FALSE, rotate = TRUE)
dev.off()
