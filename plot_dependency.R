#!/bin/Rscript

# Load the required package
library("depmap")
library("dplyr")
library("ggplot2")
library("viridis")
library("tibble")
library("gridExtra")
library("stringr")
library("depmap")
library("ExperimentHub")

# Show the package release
depmap_release()

## create ExperimentHub query object
eh <- ExperimentHub()
q <- query(eh, "depmap")
#q$preparerclass
df <- data.frame(mcols(q))
write.table(df, "entries.tsv", col.names=T, row.names=F, sep="\t")

# Get the data from 2019 Q3
#rnai <- eh[["EH3080"]]
crispr <- eh[["EH3081"]]
tpm <- eh[["EH3084"]]
metadata <- eh[["EH3086"]]

colnames(crispr)
colnames(tpm)
colnames(metadata)

# Create a histogram with the dependency scores for IKZF1
pdf("Dependency_histogram_CRISPR_Ikaros.pdf")
crispr %>% dplyr::select(gene, gene_name, dependency) %>% 
         dplyr::filter(gene_name == "IKZF1") %>% 
         ggplot(aes(x = dependency)) +
         geom_histogram(fill = '#2c7bb6') +
         geom_vline(xintercept = median(crispr$dependency, na.rm = TRUE),
                    linetype = "dotted", color = "red") +
         xlab("CERES Dependency Score (IKZF1)") + ylab("Number of DepMap cell lines")
dev.off()

# Get CRISPR dependencies per disease
meta_crispr <- metadata %>%
             dplyr::select(depmap_id, primary_disease) %>%
             dplyr::full_join(crispr, by = "depmap_id") %>%
             dplyr::filter(gene_name == "IKZF1")

# Write the CRISPR data to a file
crispr_df <- data.frame(meta_crispr)
write.table(crispr_df, "IKZF1_dependency_by_CRISPR.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Read the file
meta_crispr <- read.table("IKZF1_dependency_by_CRISPR.tsv", header=TRUE, sep="\t")
meta_crispr$primary_disease <- with(meta_crispr, reorder(primary_disease, -dependency, median, na.rm=TRUE))

# Isolate data of interest
select_line <- subset(meta_crispr, meta_crispr$cell_line=="MM1S_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
select_line

# Plot dependency across various types
pdf("IKZF1_dependency_across_various_cancer_types_crispr.pdf", useDingbats=FALSE)
p2 <- meta_crispr %>%
      ggplot(aes(x = dependency, y = primary_disease)) +
      geom_point(alpha = 0.4, size = 3.0) +
      geom_point(data=select_line, aes(x = dependency, y = primary_disease), color = 'red', size = 3.0) +
      geom_vline(xintercept=median(meta_crispr$dependency, na.rm = TRUE),
                 linetype = "longdash", color = "red") +
      xlab("CERES Dependency Score (IKZF1)") + ylab("Primary disease") + xlim(-1, 1) 
plot(p2)
dev.off()

# Get the expression for Ikaros
expression <- metadata %>%
      dplyr::select(depmap_id, primary_disease) %>%
      dplyr::full_join(tpm, by = "depmap_id") %>%
      dplyr::filter(gene_name == "IKZF1")

# Convert tribble to dataframe
expression_df <- data.frame(expression)

# Write the data to an external file 
write.table(expression_df, "IKZF1_expression_across_primary_cancer_types.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# Read the data and sort based on median expression
expression_df <- read.table("IKZF1_expression_across_primary_cancer_types.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
expression_df$primary_disease <- with(expression_df, reorder(primary_disease, expression, median, na.rm=T))

#Plot the result
pdf("IKZF1_expression_across_various_cancer_types.pdf", useDingbats=FALSE)
p2 <- ggplot(data=expression_df, aes(x = primary_disease, y = expression)) +
      geom_boxplot(outlier.alpha = 0.1, fill='#2c7bb6') +
      geom_hline(yintercept=median(expression_df$expression, na.rm = TRUE),
                 linetype = "longdash", color = "red") +
      xlab("Primary disease") + ylab("IKZF1 expression - log2(TPM+1)") +
      theme(legend.position = "none") + coord_flip()
plot(p2)
dev.off()
