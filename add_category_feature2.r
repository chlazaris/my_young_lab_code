#!/bin/Rscript

# Check for dplyr
if (!require(dplyr)) install.packages(dplyr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Check for number of inputs
if (length(args)!=3) {print("Please provide the right number of arguments. USAGE: add_category_feature2 INPUT_FILE CATEGORY FEATURE"); q(save="no")}

# Read in the file
data <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t")
colnames(data) <- c("chrom","start","end","gene","charge")

# Add category
data$category <- c(sprintf("%s", args[2]))
# Add feature
data$feature <- c(sprintf("%s", args[3]))

# Specify the name of the output
outfile <- gsub("_total.txt", "_total_charge.txt", args[1])
write.table(data, outfile, row.names=F, col.names=T, sep="\t", quote=FALSE)
