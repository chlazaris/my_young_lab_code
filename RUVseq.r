#!/bin/Rscript

# Load the required packages
library(RUVSeq)
library(RColorBrewer)

# Read in the genes
genes <- read.table("ordered_merged_gene_counts_f.tsv", header=T, sep="\t", row.names=1, check.names=F)
dim(genes)

# Filter out the genes that do not have more than 5 reads in at least two samples for each gene
filter <- apply(genes, 1, function (x) length(x[x>5])>=2)
filtered <- genes[filter,]
dim(filtered)

# Get the names of the spikes and the genes
spikes <- rownames(filtered)[grepl("^ERCC-", rownames(filtered))]
genes <- rownames(filtered)[!grepl("^ERCC-", rownames(filtered))]

# Find the number of genes and spikes that we are left with
length(genes)
length(spikes)

x <- as.factor(rep(c("0","5","10","30","240"), each=2))
y <- as.factor(rep(c(1,2), 5))
phen_data <- data.frame(row.names=colnames(filtered), Replicate=y, Time=x)
set <- newSeqExpressionSet(as.matrix(filtered), phenoData=phen_data)
pData(set)

# Inspect the data without normalization
colors <- brewer.pal(3, "Set2")

# Plot relative log expression
#pdf("RLE_plot_prenorm.pdf")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#dev.off()

# Plot PCA
#pdf("PCA_plot_prenorm.pdf")
#plotPCA(set, col=colors[x], cex=1.2)
#dev.off()

# Upper quantile normalization
set <- betweenLaneNormalization(set, which="upper")

# Plot relative log expression after normalization
pdf("RLE_plot_upperquant_norm.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

# Plot PCA after normalization
pdf("PCA_plot_upperquant_norm.pdf")
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

# Now perform normalization based on the spike-in controls
#set1 <- RUVg(set, spikes, k=1)
#pData(set1)

# Plot relative log expression after normalization
#pdf("RLE_plot_norm.pdf")
#plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#dev.off()

# Plot PCA after normalization
#pdf("PCA_plot_norm.pdf")
#plotPCA(set1, col=colors[x], cex=1.2)
#dev.off()

# Consider making groups
#differences <- makeGroups(x)
#set3 <- RUVs(set, genes, k=1, differences)
#pData(set3)

# Plot relative log expression after normalization with groups
#pdf("RLE_plot_group_norm.pdf")
#plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#dev.off()

# Plot PCA after normalization with groups
#pdf("PCA_plot_group_norm.pdf")
#plotPCA(set3, col=colors[x], cex=1.2)
#dev.off()

# Differential expression analysis
# How to set up the model
# https://support.bioconductor.org/p/113630/

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set), colData = pData(set), design = ~ Replicate + Time)
dds <- DESeq(dds, test="LRT", reduced = ~ Replicate)
time_5_vs_0 <- data.frame(results(dds, contrast=c("Time", "5", "0")))
time_10_vs_0 <- data.frame(results(dds, contrast=c("Time", "10", "0")))
time_30_vs_0 <- data.frame(results(dds, contrast=c("Time", "30", "0")))
time_240_vs_0 <- data.frame(results(dds, contrast=c("Time", "240", "0")))

write.table(time_5_vs_0, "expr_5_vs_0.tsv", row.names=T, col.names=T, sep="\t", quote=F)
write.table(time_10_vs_0, "expr_10_vs_0.tsv", row.names=T, col.names=T, sep="\t", quote=F)
write.table(time_30_vs_0, "expr_30_vs_0.tsv", row.names=T, col.names=T, sep="\t", quote=F)
write.table(time_240_vs_0, "expr_240_vs_0.tsv", row.names=T, col.names=T, sep="\t", quote=F)

#dexpr_data <- data.frame(res)
#write.table(expr_data, "diff_expr_data.tsv", row.names=T, col.names=T, sep="\t", quote=F)

#dds <- DESeq(dds, test="LRT", reduced=as.formula("~ W_1"))
#res <- results(dds)

# Export the data
#expr_data <- data.frame(res)
#colnames(res) <- resultsNames(dds)
#write.table(expr_data, "diff_expr_data.tsv", row.names=T, col.names=T, sep="\t", quote=F)
