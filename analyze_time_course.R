# Use the required libraries
library(DESeq2)

# Check for the number of arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
    print("Please provide right number of arguments...")
    print("USAGE: run_deseq2 COUNT_MATRIX SAMPLE_SHEET")
    q(save="no")
}

# Read in the table 
data <- read.table(sprintf("%s", args[1]), header=T, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

# Reorder the data based on what is in the sample sheet
col_order <- c(1,4,11,10,8,9,3,7,5,2,6)
data <- data[, col_order]
head(data)

# Convert the data to count matrix
countdata <- as.matrix(data[2:ncol(data)])
rownames(countdata) <- data$Gene
head(countdata)

# Read in the sample sheet
coldata <- read.table(sprintf("%s", args[2]), header=T)
coldata$rep <- as.factor(coldata$rep)
coldata$time <- as.factor(coldata$time)

# Run DEseq2
ddsMat<-DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design = ~rep + time)

# Exclude the ones that do not have more than 1 count (e.g. all zeros or the ones with only 1 in one of the samples)
ddsMat <- ddsMat[ rowSums(counts(ddsMat)) > 1, ]

# Find the ones that are differentially expressed compared to timepoint 0
ddsLRT <- DESeq(ddsMat, test="LRT", reduced = ~rep)

# Get the results for all the comparisons
comparisons <- resultsNames(ddsLRT)
comparisons <- comparisons[3:length(comparisons)]
comparisons

# Save all the significant results for the individual comparisons (all time points) 
for (comparison in comparisons) {
	res <- results(ddsLRT, name=sprintf("%s", comparison), test="Wald")
	res$Gene <- row.names(res)
	sig_res <- subset(res, res$padj < 0.05)
	outfile=paste0(comparison,"_diff_expressed.tsv")
	write.table(sig_res, outfile, row.names=F, col.names=T, sep="\t", quote=F)
}
