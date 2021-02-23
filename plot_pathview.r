#!/bin/Rscript

# Download PathView
#if (!requireNamespace("BiocManager", quietly=TRUE))
#install.packages("BiocManager")
#BiocManager::install("pathview")

# Load PathView
library(pathview)
library(data.table)
library(dplyr)

# Check the arguments
args <- commandArgs(trailingOnly=TRUE)

# Read in the gene symbols and the pathways
gene_symbols <- read.table(sprintf("%s", args[1]), header=F, sep="\t", stringsAsFactors=FALSE)
kegg_pathways <- fread(sprintf("%s", args[2]), header=F, sep="\t", colClasses = 'character', data.table = FALSE)

# Give headers
colnames(gene_symbols) <- c("symbol")
colnames(kegg_pathways) <- c("pathway")

pathways <- as.vector(as.character(kegg_pathways$pathway)) 
symbols <- as.vector(as.character(gene_symbols$symbol))

# Convert the pathway names to compare with hsa pathways
paths <- paste0("hsa", pathways)

# Isolate only the pathways that are present in KEGG
data(paths.hsa)
pathways <- gsub("hsa", "", paths[paths %in% names(paths.hsa)])

# Conversion of gene symbol to EntrezID
entrez_ids <- id2eg(symbols, org = "Hs", pkg.name = "org.Hs.eg.db")

# Create a vector with 1s and names 
# the Entrez IDs
gene_ids <- rep(1, length(symbols))
names(gene_ids) <- entrez_ids[,2]

# Create a list to store all the data for the SE-regulated genes (percentage) for all pathways
l <- list()

# Check if the pathways of interest are in the whole database of pathways
for (i in 1:length(pathways)) {
	pv.out <- pathview(gene.data = gene_ids, gene.idtype="entrez", 
		   	pathway.id =  pathways[i], 
		   	species = "hsa", 
                   	kegg.native = T)
	# Calculate the fraction of super-enhancer regulated genes and save it with the pathway name
	df <- pv.out$plot.data.gene 
	df$color_value <- ifelse(df$mol.col=="#FF0000", 1 , 0)
	print(df)

	# Find unique genes involved in the pathway and how many are super-enhancer regulated (isSuper = 1)
	data <- df %>% select(kegg.names, labels, color_value) %>% group_by(kegg.names, labels) %>% summarize(isSuper=max(color_value)) %>% data.frame()
	#print(data)
	data$pathway <- rep(paste0("hsa", pathways[i]), nrow(data))
	colnames(data) <- c("kegg_name", "geneSymbol", "isSuper", "pathway")
	print(data)

	# Save the data into a file
	outname <- paste0("hsa", pathways[i], "_genes.txt")
	write.table(data, outname, row.names=F, col.names=T, sep="\t", quote=F)

	# Get the super-enhancer count (unique genes in red)
	se_count <- sum(data$isSuper)

	# Get the total gene count (unique genes in the pathway)
	total_count <- length(data$isSuper)
        
	# Get the percentage of super-enhancer regulated genes over the total 
	percent <- as.numeric(as.character((se_count/total_count)*100))
        pathway_name <- paste0("hsa", pathways[i])
	
	# Store all these elements to the corresponding position in the list
	l[[i]] <- c(pathway_name, se_count, total_count, percent)
}

# Convert all the list elements to a dataframe with the percentages of se_genes
df_se_percent <- data.frame(do.call("rbind", l))
colnames(df_se_percent) <- c("pathway", "se_genes", "total_genes", "percentage")
df_se_percent$percentage <- round(as.numeric(as.character(df_se_percent$percentage)), 2)
df_se_percent$se_genes <- as.numeric(as.character(df_se_percent$se_genes))
df_se_percent$total_genes <- as.numeric(as.character(df_se_percent$total_genes))
write.table(df_se_percent, "SE_percentage.tsv", col.names=T, row.names=F, sep="\t", quote=F)
