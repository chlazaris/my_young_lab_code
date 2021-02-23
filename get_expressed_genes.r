#!/bin/Rscript

# This is a script that gets the gene abundance
# output from the RNA-seq pipeline and calculates 
# gene size in kb and log2FPKM

# Load the required packages
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(dplyr)) install.packages("dplyr")

# Load the required libraries
library(ggpubr)
library(dplyr)

# Get the input
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]

# Read in the input file
data <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Get the expressed genes with FPKM greater than 1
data <- subset(data, data$FPKM > 1)

# Get the gene size in kb, log(10)FPKM and log(2)FPKM
data2 <- data %>% mutate(gene_length=(End-Start+1), gene_length_kb=((End-Start+1)/1000), logFPKM=log(FPKM), log2FPKM=log2(FPKM)) %>% data.frame()

# Get only the genes that are longer than 1Kb
data2 <- subset(data2, gene_length > 1000)

# Get the output file
outfile <- gsub(".txt", "_expressed.tsv", args[1])

# Write the data to a file
write.table(data2, outfile, row.names=F, col.names=T, sep="\t", quote=F)

