#!/bin/Rscript

# Load the required libraries
library("LOLA")
library("simpleCache")
library("GenomicRanges")

regionDB = loadRegionDB("/lab/solexa_young/lazaris/data/LOLA/LOLACore/hg19")
areas_of_interest <- readBed("SEs.bed")
activeDHS <- readBed("/lab/solexa_young/lazaris/data/LOLA/lola_vignette_data/activeDHS_universe.bed")

# Convert to GRanges object
userSets <- GRangesList(areas_of_interest)

# Get the results of the Fisher's exact test (enrichment)
locResults <- runLOLA(userSets, activeDHS, regionDB, cores=1)

# Export the data into a file
write.table(locResults,"LOLA_results.tsv",col.names=T,row.names=F,sep="\t",quote=F)
