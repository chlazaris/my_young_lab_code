#!/bin/Rscript 

# Load the required packages
# if (!require("fabricatr")) {install.packages("fabricatr")}
library(dplyr)
library(ggpubr)

# Load the required libraries
# library(fabricatr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Read in the data
data <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Exclude the genes that are not expressed
exprs <- subset(data, data$FPKM > 0)
exprs$gene_size_kb <- ((exprs$End - exprs$Start) + 1)/1000

# Get the quantiles
quants <- mutate(exprs, quantile_rank = ntile(exprs$FPKM, 5), size_quantile_rank = ntile(exprs$gene_size_kb, 5)) 

# Set the quantiles as factors
quants$quantile_rank <- as.factor(quants$quantile_rank)
quants$size_quantile_rank <- as.factor(quants$size_quantile_rank)
head(quants)
str(quants)

pdf("test.pdf")
ggboxplot(quants, x = "size_quantile_rank", y = "FPKM",
          color = "size_quantile_rank", palette = "jco")+ yscale("log10") + stat_compare_means()
dev.off()
