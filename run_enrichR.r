#!/usr/bin/env Rscript

# Load the required libraries
library(enrichR)
setEnrichrSite("Enrichr") # Human genes

# List the available databases
dbs <- listEnrichrDbs()
dbs2 <- dbs
dbs2

# Get only the ones that are online
websiteLive <- TRUE
if (is.null(dbs)) websiteLive <- FALSE

# Select your favorite databases
dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", 
	 "GO_Molecular_Function_2021", "WikiPathway_2021_Human", 
	 "KEGG_2021_Human", "HDSigDB_Human_2021", 
	 "MSigDB_Hallmark_2020", "DSigDB", "Reactome_2016", "Elsevier_Pathway_Collection")

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
	print("Please provide the correct number of arguments...")
	print("run_enrichR.r [GENE LIST; .txt]")
}

# Get the list of genes
genes <- read.table(sprintf("%s", args[1]), header=T, sep="\t", check.names=F)
colnames(genes) <- c("Genes")
genes <- as.vector(as.character(genes$Genes))

# Select the genes that you want to check
if (websiteLive) {
    enriched <- enrichr(genes, dbs)
}

# Go through all the databases
# and get the significant hits
for (i in 1:length(dbs)) {
	db_name <- dbs[i]
	result <- enriched[[i]]
	
	# Write the significant hits in a file
	outfile <- paste0(db_name, "_hits.tsv")
	write.table(result, outfile, row.names=F, col.names=T, quote=F, sep="\t")
	
	# Plot the top 20 ones from each database
	outfile1 <- paste0(db_name, "_top20_hits.pdf")
	pdf(outfile1)
	p <- plotEnrich(result, showTerms = 20, numChar = 60, y = "Ratio", orderBy = "Adjusted.P.Value")
	print(p)
	dev.off()
}
