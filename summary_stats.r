#!/bin/Rscript

# Install the required packages
if (!require("dplyr")) install.packages("dplyr")

# Load the required packages
library(dplyr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Check if the number of arguments is correct
if (length(args)!=1) {print("Please provide right number of arguments"); print("USAGE: summary_stats INPUT_FILE"); q(save="no")}

# Store the name for the input
input <- sprintf(args[1])

# Read in the data
data <- read.table(sprintf("%s", args[1]), header=T, sep="\t", stringsAsFactors=FALSE)

# Make sure to get the entries where RPK > 0
data <- subset(data, data$RPK > 0)

# Summarize the data
summary <- data %>% group_by(category, feature) %>% summarise(n=n(), min=min(RPK), mean=mean(RPK), median=median(RPK), max=max(RPK), std=sd(RPK)) %>% data.frame()
summary

# Store to output
outfile <- gsub(".tsv", "_summary_stats.tsv", input)
write.table(summary, outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

