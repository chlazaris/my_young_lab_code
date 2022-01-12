#!/bin/env Rscript

# This is a script that takes a peak file and the annotations for the peaks
# and plots the signal, the quantiles and where peaks for certain genes fall
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

# Get the file with the signal for peaks
signal <- read.table(args[1], header=F, sep="\t", check.names=F)
head(signal)
colnames(signal) <- c("chrom","chromStart","chromEnd","PeakID",
		      "score","strand","signalValue","pValue",
		      "qValue","peak")

# Get the gene annotations
annots <- data.frame(read_tsv(args[2]))
colnames(annots)[1] <- c("PeakID")
head(annots)

# Join the two
joined_data <- inner_join(signal, annots, by="PeakID")
joined_data$rank <- c(1:nrow(joined_data))
write.table(joined_data,"annotated_peak_data_with_signal_rank.tsv", sep="\t", quote=F, row.names=F, col.names=T)

# Get the histogram
# Get the quantiles and plot them
quantiles=quantile(joined_data$signalValue)
sink("quantiles.txt")
print(quantiles)
sink()

pdf("signal_histogram_with_quantiles.pdf")
hist(joined_data$signalValue, n=100, xlab="Peak signal", ylab="Frequency")
abline(v = quantiles, col='red', lty=3)
dev.off()

# Plot lines for certain genes
select_data <- subset(joined_data, joined_data$Gene.Name %in% c("TIMM22", "LARS", "GSK3A", "FTO", "PANK3", "TSC2", "CAPN10"))
select_signal <- select_data$signalValue

pdf("signal_histogram_with_quantiles_and_certain_genes_highlighted.pdf")
hist(joined_data$signalValue, n=100, xlab="Peak signal", ylab="Frequency")
abline(v = quantiles, col='red', lty=3)
abline(v = select_signal, col='blue', lty=3)
dev.off()

# Plot lines for genes including FASN, SMAD7
select_data2 <- subset(joined_data, joined_data$Gene.Name %in% c("FASN", "SMAD7"))
select_signal2 <- select_data2$signalValue

pdf("signal_histogram_with_quantiles_and_genes_highlighted.pdf")
hist(joined_data$signalValue, n=100, xlab="Peak signal", ylab="Frequency")
abline(v = quantiles, col='black', lty=2)
abline(v = select_signal, col='blue', lty=3)
abline(v = select_signal2, col='red', lty=3)
dev.off()
