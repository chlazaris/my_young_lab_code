#!/bin/Rscript

# Load the required libraries
library("ggpubr")
library("dplyr")

# Get the input .csv file
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("Provide the right number of arguments..."); q(save="no")}
input <- args[1]

# Read in the data
data <- read.table(sprintf("%s", input), header=F, sep="\t")
colnames(data) <- c("gene", "count", "total")

# Calculate percentage
data <- data %>% mutate(percent=round(count/total * 100, 1)) %>% data.frame()

# Get the name of output
output <- gsub("tsv","pdf",input)

# Print the plot
pdf(output)
p <- ggbarplot(data, "gene", "percent",
          fill = "blue",               # fill color
          sort.val = "desc",           # Sort the value in descending order
          sort.by.groups = FALSE,      # Don't sort inside each group
          x.text.angle = 45,           # Rotate vertically x axis texts
	  label = FALSE, lab.pos = "in", lab.col = "white", lab.size=1.5
	  )
p + font("xy.text", size = 10) + scale_y_continuous(expand = c(0, 0)) + xlab("Gene") + ylab("Frequency (%)")
dev.off()
