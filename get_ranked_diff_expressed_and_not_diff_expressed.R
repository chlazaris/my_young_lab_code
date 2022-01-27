#!/bin/Rscript

# Load the required libraries
library(readr)

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
	print("Please provide the correct argument number...")
	print("USAGE: get_ranked_diff_expressed_and_not_diff_expressed.R [DEseq2 result]")
	q(save="no")
}

data <- read.table("time_240_vs_0_deseq2_results.tsv")
data$gene <- row.names(data)
dim(data)

# Get differentially expressed genes
diff_expr <- subset(data, data$padj < 0.05)
order_diff_expr <- diff_expr[order(-(abs(diff_expr$log2FoldChange))),]
#head(order_diff_expr)
#tail(order_diff_expr)

# Get the genes that are not differentially expressed
not_diff_expr <- subset(data, data$padj > 0.05)
order_not_diff_expr <- not_diff_expr[order(-(abs(not_diff_expr$log2FoldChange))),]
#head(order_not_diff_expr)
#tail(order_not_diff_expr)

# Get both in the same data frame
t <- rbind(order_diff_expr, order_not_diff_expr)
write.table(t,"diff_and_not_diff_expressed_ranked_by_log2FC.tsv", row.names=F, col.names=T, sep="\t", quote=F)
