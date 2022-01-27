#!/bin/Rscript

# Load the required libraries
#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz', repos=NULL)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("affycoretools")

#library(ff)
library(oligo)
library(GEOquery)
library(affycoretools)
library(limma)
library(dplyr)

## ---getData---
#getGEOSuppFiles("GSE31365")
#untar("GSE31365/GSE31365_RAW.tar", exdir = "GSE31365/CEL")

## ---Read the data---
celfiles <- list.files("GSE31365/CEL", full=TRUE)
rawData <- read.celfiles(celfiles)

## --phenoData--
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub("^[^_]*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("JQ1", sampleNames(rawData)), "JQ1", "DMSO")
pData(rawData)

## --rma--
normData <- oligo::rma(rawData, normalize=T, background=T)
#normData

# Annotate Data
normData <- annotateEset(normData, columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"), type="core", pd.hugene.1.0.st.v1)

## --Plot the norm data--
pdf("norm_data_boxplot.pdf")
boxplot(exprs(normData))
dev.off()

# Subset to focus only on MM1S
mmData <- normData[,sampleNames(rawData) %in% c("MM1S_DMSO_A","MM1S_DMSO_B","MM1S_JQ1_A","MM1S_JQ1_B")]
mmData$treatment <- factor(c("DMSO","DMSO","JQ1","JQ1"))

# Apply the limma model
design <- model.matrix(~mmData$treatment - 1)
colnames(design) <- c("DMSO", "JQ1")

# Use limma (linear modeling)
fit <- lmFit(mmData, design)
contrast.matrix <- makeContrasts("JQ1-DMSO", levels=design)
contrast.matrix
fitC <- contrasts.fit(fit, contrast.matrix)
fitC <- eBayes(fitC)

# Get the differentially expressed genes using lfc=1 and adjusted p-value<0.05
de_genes <- topTable(fitC, coef=1, lfc=1, number=Inf, p.value=0.05, adjust.method="BH")
# Exclude the probe ID
de_genes <- de_genes[,2:ncol(de_genes)]
# Remove the duplicated raws
de_genes <- de_genes[!duplicated(de_genes),]
# Romove the ones where the identifier in not available
de_genes <- data.frame(subset(de_genes, !is.na(de_genes$ID)))

# Select the one with the highest average expression when the gene is the same
de_genes <- de_genes %>% group_by(SYMBOL) %>% top_n(n=1, wt=AveExpr) %>% data.frame()

# Get the top 50 upregulated and top 50 downregulated genes
top_genes <- topTable(fitC, coef=1, number=108, adjust.method="BH")
#top_genes <- top_genes[order(-top_genes$logFC),]
#top_genes <- data.frame(subset(top_genes, !is.na(top_genes$ID)))
#top_genes <- top_genes %>% group_by(SYMBOL) %>% top_n(n=1, wt=AveExpr) %>% data.frame()
#head(top_genes, n=50)
#tail(top_genes, n=50)
#top_100 <- rbind(head(top_genes, n=50), tail(top_genes, n=50))

# Write to a file
write.table(de_genes, "diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)
write.table(top_genes, "top_100_diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)
