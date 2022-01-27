#!/bin/Rscript

# This is a script that runs DiffBind
# accepting as input the corresponding sample sheet

# Load the required packages and install
# if they are not present
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DiffBind)

# Check the number of arguments
input <- commandArgs(trailingOnly=TRUE)

if (length(input)!=1) {
	print("Please provide the right number of arguments...")
	print("USAGE: Rscript run_diffbind.r [SAMPLE_SHEET (.csv)]")
	q(save="no")
}

# Read in the sample sheet
data <- read.csv(sprintf("%s", input[1]), header=TRUE, stringsAsFactors=FALSE)
samples <- dba(sampleSheet=data)

# Get the Factor of interest
factor <- unique(data$Factor)

# Create a correlation map
# to show the relationship of the samples
outname <- paste(factor, "correlation_plot.pdf", sep="_")
pdf(outname, useDingbats=FALSE)
plot(samples)
dev.off()

# Count reads and FRiP score for each one of the samples
reads <- dba.count(samples)

# Normalize the reads based on library sizes
#norm <- dba.normalize(samples, bRetrieve=TRUE)
norm <- reads

# Create the contrasts
contrasts <- dba.contrast(norm, categories=DBA_CONDITION, minMembers=2)
contrasts <- dba.analyze(contrasts)
dba.show(contrasts, bContrasts=TRUE)

# Plot a heatmap showing only the differentially bound sites
outname1 <- paste(factor, "correlation_heatmap_based_on_differentially_bound_sites.pdf", sep="_")
pdf(outname1, useDingbats=FALSE)
plot(contrasts, contrast=1)
dev.off()

# Retrieve the differentially bound sites
outname2 <- paste(factor, "db_sites.tsv", sep="_")
db_sites <- dba.report(contrasts)
df <- data.frame(db_sites)
write.table(df, outname2, row.names=F, col.names=T, sep="\t", quote=F)

# PCA plot
outname3 <- paste(factor, "pca_plot.pdf", sep="_")
pdf(outname3, useDingbats=FALSE)
dba.plotPCA(contrasts, attributes=c(DBA_CONDITION), label=DBA_TISSUE)
dev.off()

# Get the results with edgeR and DESeq2
results <- dba.analyze(contrasts, method=DBA_ALL_METHODS)
dba.show(results, bContrasts=TRUE)

# Get a plot of the overlap between edgeR and DESeq2 calls
outname4 <- paste(factor, "overlap_edgeR_DESeq2_plot.pdf", sep="_")
pdf(outname4, useDingbats=FALSE)
dba.plotVenn(results, contrast=1, method=DBA_ALL_METHODS, bDB=TRUE)
dev.off()
