#!/bin/Rscript
library(dplyr)
library(reshape2)

# Import the data
data <- read.table("select_SE_genes_per_patient.tsv", header=F, sep="\t")

# Name the columns
colnames(data) <- c("donor", "se_gene")

# Add a column with the counts 
data2 <- data %>% group_by(donor, se_gene) %>% mutate(count=n()) %>% as.data.frame()

# Create the matrix to plot
df <- dcast(data2, donor~se_gene)
df[is.na(df)] <- 0

# Convert to matrix
m <- as.matrix(df[,2:ncol(df)])
rownames(m) <- df$donor

# Specify the colors
color_palette <- c('#ef8a62', '#67a9cf')

# Create the heatmap
pdf("heatmap.pdf", height=25, width=10, useDingbats=FALSE)
heatmap(m, scale="none", col=rev(color_palette))
dev.off()



