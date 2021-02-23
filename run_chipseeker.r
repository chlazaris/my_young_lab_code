#!/bin/Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("stephenturner/annotables")

# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
#library(annotables)
library(org.Hs.eg.db)

# Load data
samplefiles <- list.files(".", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("IKZF1")

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 1000), verbose=FALSE)

# Plot the percentages of binding in features
#pdf("binding_profile_piechart.pdf")
#plotAnnoPie(peakAnnoList)
#dev.off()

pdf("binding_profile_barplot.pdf")
plotAnnoBar(peakAnnoList)
dev.off()

pdf("binding_profile_distance_from_TSS")
plotDistToTSS(peakAnnoList)
dev.off()

#peak_annot_percent <- as.data.frame(peakAnnoList[[1]])
