#!/bin/Rscript

# This is to create a grouped violin and
# also plot the stats

# Install and load the required libraries
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)

# Check the arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("Please provide correct number of arguments"); q(save="no")}
input <- args[1]

# Read the input data
data <- read.table(sprintf("%s", input), header=T, sep="\t", stringsAsFactors=FALSE)

# Convert the correct columns to factors
data$category <- as.factor(data$category)
data$feature <- factor(data$feature, levels=c("enhancer", "proximal", "gene"))

# Get only the data where RPK>0
data <- subset(data, data$RPK > 0)

# Create the violin plot
pdf("test.pdf", useDingbats = FALSE)
p <- ggviolin(data, "feature", "RPK", fill = "category",
palette = c("#D3D3D3", "#ffffff"), add = "median") + yscale("log10", .format = TRUE) + ylab("Reads per kb (RPK)") + xlab("Feature") 
p + stat_compare_means(aes(group = category))
print(p)
dev.off()

#stat_compare_means(aes(group = category), data = data)
#my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
#p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#  stat_compare_means(label.y = 50)                   # Add global p-value

