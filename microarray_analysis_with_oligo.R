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
normData <- annotateEset(normData, type="core", pd.hugene.1.0.st.v1)

## --Plot the norm data--
pdf("norm_data_boxplot.pdf")
boxplot(exprs(normData))
dev.off()

# Subset to focus only on MM1S
mmData <- normData[,sampleNames(rawData) %in% c("MM1S_DMSO_A","MM1S_DMSO_B","MM1S_JQ1_A","MM1S_JQ1_B")]
mmData$treatment <- factor(c("DMSO","DMSO","JQ1","JQ1"))

# Apply the limma model
design <- model.matrix(~mmData$treatment - 1)
design 
colnames(design) <- c("DMSO", "JQ1")
fit <- lmFit(mmData, design)
contrast.matrix <- makeContrasts("JQ1-DMSO", levels=design)
contrast.matrix
fitC <- contrasts.fit(fit, contrast.matrix)
fitC <- eBayes(fitC)
de_genes <- topTable(fitC, coef=1, lfc=1, number=Inf, p.value=0.05, adjust.method="BH")
de_genes <- de_genes[,2:ncol(de_genes)]
de_genes <- de_genes[!duplicated(de_genes),]
write.table(de_genes, "diff_expressed_genes.tsv", row.names=F, col.names=T, sep="\t", quote=F)

## --write to file
#df <- exprs(normData)[1:nrow(normData),1:ncol(normData)]
#head(df)
#pdata <- pData(eset)
#d <- cbind(pdata, t(m))
