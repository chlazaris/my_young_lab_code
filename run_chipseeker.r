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
samplefiles <- list.files("peaks/", pattern= ".bed", full.names=T)
#samplefiles <- as.list(samplefiles)
samplefiles <- samplefiles[1]
samplefiles
#names(samplefiles) <- c("peaks")
type(samplefiles)

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get annotations (define promoters as TSS +/-2kb)
peakAnno <- annotatePeak(samplefiles, tssRegion=c(-2000, 2000), TxDb=txdb, verbose=FALSE, annoDb="org.Hs.eg.db")
annotations <- as.data.frame(peakAnno)
write.table(annotations, "peak_annotations.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# Plot the percentages of binding in features
pdf("binding_profile_piechart.pdf")
plotAnnoPie(peakAnno)
dev.off()

pdf("binding_profile_barplot.pdf")
plotAnnoBar(peakAnno)
dev.off()

pdf("binding_profile_vennpie.pdf")
vennpie(peakAnno)
dev.off()

pdf("binding_profile_vennpie_upsetplot.pdf")
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

pdf("binding_profile_distance_from_TSS.pdf")
plotDistToTSS(peakAnno)
dev.off()
