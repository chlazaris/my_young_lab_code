#!/bin/Rscript

# Install BioMart
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load the required libraries
library(biomaRt)

# Use the right version from Ensembl
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
ensembl67=useMart(host="may2012.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl",ensembl67)

# Get the appropriate attributes (Ensembl ID, Entrez ID, Gene Symbol)
# mapping <- getBM(attributes=c('refseq_mrna', 'ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), mart=mart)
mapping <- getBM(attributes=c('refseq_mrna', 'ensembl_gene_id', 'external_gene_id'), mart=mart)
refseq <- mapping[!is.na(mapping$refseq_mrna),]
head(refseq)
write.table(refseq, "gene_mapping_mm9.tsv", row.names=F, col.names=T, sep="\t", quote=F)
