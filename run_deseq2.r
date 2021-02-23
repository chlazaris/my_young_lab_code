#!/bin/Rscript

# This is a script to run DESeq2 using featureCounts
# matrix as input

###############################################################################################################
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
################################################################################################################

# Use the required libraries
library(DESeq2)

# Check for the number of arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
	print("Please provide right number of arguments...")
    print("USAGE: run_deseq2 COUNT_MATRIX SAMPLE_SHEET"); q(save="no")
}

# Read in the table 
data <- read.table(sprintf("%s", args[1]), header=T, sep="\t", check.names=F, row.names=1)

# Convert the data to count matrix
countdata <- as.matrix(data)

# Read in the sample sheet
coldata <- read.table(sprintf("%s", args[2]), header=T, sep="\t", check.names=F, row.names=1)

# Convert the condition column to factor
coldata$condition <- factor(coldata$condition, levels=c("DMF","cisplatin"))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)

# Get the results
res <- results(dds)
head(res)

# Combine with normalized data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- c("Gene")
# Get all the results
write.table(resdata, "all_deseq2_results.tsv", col.names=TRUE, row.names=FALSE, quote=F, sep="\t")

# Get differentialy expressed genes
diff_exp <- subset(resdata, abs(resdata$log2FoldChange)>1 & resdata$padj<0.05)
head(diff_exp)
write.table(diff_exp, "diff_expr_log2FC_1_padj_005.tsv", col.names=T, row.names=T, quote=F, sep="\t")

# Remove the ones with NA values
# res2 <- subset(res, !is.na(res$log2FoldChange))

# Create a volcano plot
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=30)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-3, 3))
dev.off()





