#!/bin/Rscript

# Perform GO analysis for cis-regulatory elements
library(rGREAT)

# Install it if it is not there 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#	install.packages("BiocManager")
#BiocManager::install("rGREAT")

# Read in the arguments
args <- commandArgs(trailingOnly=TRUE)

# Read the .bed file and the reference genome (hg19, mm9)
bed <- read.table(sprintf("%s", args[1]), header=F, sep="\t")
colnames(bed) <- c("chr", "start", "end", "description", "score", "strand")
organism <- as.vector(as.character(args[2]))

# Run the job using the latest version of the tool
job = submitGreatJob(bed, species = organism, version = "4.0.4")

# Get the enrichment tables for all available ontologies
all_ontologies <- availableOntologies(job)
tb = getEnrichmentTables(job, ontology = all_ontologies)

# Write the outputs in separate files
for (i in 1:length(all_ontologies)) {
	name <- gsub(" ", "", all_ontologies[i])
	data <- tb[[i]]
	write.table(data, name, col.names=T, row.names=F, sep="\t", quote=F)
}

