#!/bin/Rscript

# Install DepMap package
#remotes::install_github("UCLouvain-CBIO/depmap", ref = "4a9f52ed6cf9c3821891ebdd9db317194d6518c9")

# Load the required package
library("dplyr")
library("ggplot2")
library("viridis")
library("tibble")
library("gridExtra")
library("stringr")
library("depmap")
library("ExperimentHub")

# Check the number of arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=1 ){
	print("Please provide the right number of arguments...")
	print("USAGE: plot_dependency.R [GENE NAME]")
	q(save="no")
}

# Get the gene name
genename <- as.character(args[1])

# Show the package release
depmap_release()

## create ExperimentHub query object
eh <- ExperimentHub()
q <- query(eh, "depmap")
df <- data.frame(mcols(q))
write.table(df, "entries.tsv", col.names=T, row.names=F, sep="\t")

# Get the data from 2019 Q3
#rnai <- eh[["EH3080"]]
crispr <- eh[["EH3081"]]
tpm <- eh[["EH3084"]]
metadata <- eh[["EH3086"]]

# Specify the name of the label 
x_label <- sprintf("CERES Dependency Score (%s)", genename)

# Create a histogram with the dependency scores for IKZF1
pdf(sprintf("Dependency_histogram_CRISPR_%s.pdf", genename))
crispr %>% dplyr::select(gene, gene_name, dependency) %>% 
         dplyr::filter(gene_name == sprintf("%s", genename)) %>% 
         ggplot(aes(x = dependency)) +
         geom_histogram(fill = '#2c7bb6') +
         geom_vline(xintercept = median(crispr$dependency, na.rm = TRUE),
                    linetype = "dotted", color = "red") +
         xlab(x_label) + ylab("Number of DepMap cell lines")
dev.off()

# Get CRISPR dependencies per disease
meta_crispr <- metadata %>%
             dplyr::select(depmap_id, primary_disease) %>%
             dplyr::full_join(crispr, by = "depmap_id") %>%
             dplyr::filter(gene_name == sprintf("%s", genename))

# Write the CRISPR data to a file
crispr_df <- data.frame(meta_crispr)
outfile <- paste(genename, "dependency_by_CRISPR.tsv", sep = "_")
write.table(crispr_df, outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Read the file
meta_crispr <- read.table(sprintf("%s", outfile), header=TRUE, sep="\t")
meta_crispr$primary_disease <- with(meta_crispr, reorder(primary_disease, -dependency, median, na.rm=TRUE))

# Get the output
outfile1 <- paste(genename, "dependency_across_various_cancer_types_crispr.pdf", sep = "_")

# Plot dependency across various types
pdf(outfile1, useDingbats=FALSE)
p2 <- meta_crispr %>%
      ggplot(aes(x = dependency, y = primary_disease)) +
      geom_point(alpha = 0.4, size = 3.0) +
      geom_vline(xintercept=median(meta_crispr$dependency, na.rm = TRUE),
                 linetype = "longdash", color = "red") +
      xlab(x_label) + ylab("Primary disease") 
plot(p2)
dev.off()

# Get the expression for the gene of interest
expression <- metadata %>%
      dplyr::select(depmap_id, primary_disease) %>%
      dplyr::full_join(tpm, by = "depmap_id") %>%
      dplyr::filter(gene_name == sprintf("%s", genename))

# Convert tribble to dataframe
expression_df <- data.frame(expression)

outfile2 <- paste(genename, "expression_across_primary_cancer_types.tsv", sep="_")
# Write the data to an external file 
write.table(expression_df, outfile2, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# Read the data and sort based on median expression
expression_df <- read.table(outfile2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
expression_df$primary_disease <- with(expression_df, reorder(primary_disease, expression, median, na.rm=T))

outfile3 <- paste(genename, "expression_across_various_cancer_types.pdf", sep="_")
y_label <- sprintf("%s expression - log2(TPM+1)", genename)

#Plot the result
pdf(outfile3, useDingbats=FALSE)
p2 <- ggplot(data=expression_df, aes(x = primary_disease, y = expression)) +
      geom_boxplot(outlier.alpha = 0.1, fill='#2c7bb6') +
      geom_hline(yintercept=median(expression_df$expression, na.rm = TRUE),
                 linetype = "longdash", color = "red") +
      xlab("Primary disease") + ylab(y_label) +
      theme(legend.position = "none") + coord_flip()
plot(p2)
dev.off()
