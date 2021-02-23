#!/bin/Rscript

# Load the required libraries
library(dplyr)
library(ggpubr)

data <- read.table("all_enhancers.txt", sep="\t", header=T)

# Select the gene column only and find the gene occurrences per enhancer 
data2 <- data %>% select(gene) %>% group_by(gene) %>% mutate(count=n())
data3 <- distinct(data2) %>% data.frame()

total_genes <- dim(data3)[1]
genes_with_unique_enhancer <- dim(data3[data3$count==1,])[1]
genes_with_more_enhancers <- dim(data3[data3$count>1,])[1]

pct_unique <- round((genes_with_unique_enhancer/total_genes)*100, 2)
pct_multiple <- round((genes_with_more_enhancers /total_genes)*100, 2)

# pct_unique: 68.59
# pct_multiple: 31.41

# Plot a histogram with the frequencies
pdf("test.pdf")
gghistogram(data3, x = "count", bins = 10, 
            fill = "#0073C2FF", color = "#0073C2FF")
dev.off()
