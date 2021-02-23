#!/bin/Rscript

# Get the required packages
if (!require(tidyverse)) install.packages("tidyverse", repos='http://cran.us.r-project.org')
if (!require(dplyr)) install.packages("dplyr", repos='http://cran.us.r-project.org')

# Load the required libraries
library(ggplot2)
library(dplyr)

# Check for the input file
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {print("Please provide correct number of arguments"); q(save="no")}

# Read in the input file
data <- read.table(sprintf("%s", args[1]), header=TRUE, stringsAsFactors=FALSE)

# Create the boxplot
p <- ggplot(data, aes(x=feature, y=RPK, fill=feature)) + ylab("Reads per kb (RPK)") + geom_violin() + geom_boxplot(width=0.1, fill="white", outlier.size=0.3) + scale_y_continuous(trans='log10') + facet_wrap(~category)

# Create the output file
outfile <- gsub(".tsv", "_boxplot.pdf", args[1])

# Print
pdf(sprintf("%s", outfile))
print(p)
dev.off()

# Give summary statistics
summary <- data %>% group_by(category, feature) %>% summarise(n=n(), Min=min(RPK), Max=max(RPK), Mean=mean(RPK), Median=median(RPK), Sd=sd(RPK)) %>% data.frame()
summary
write.table(summary, "summary_stats.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)

# Give summary statistics (for no zero RPK)
data_nozero_rpk <- subset(data, data$RPK > 0)
summary_nozero_rpk <- data_nozero_rpk %>% group_by(category, feature) %>% summarise(n=n(), Min=min(RPK), Max=max(RPK), Mean=mean(RPK), Median=median(RPK), Sd=sd(RPK)) %>% data.frame()
write.table(summary, "summary_stats_nozero_rpk.tsv", col.names=T, row.names=F, sep="\t", quote=FALSE)

