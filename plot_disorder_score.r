#!/bin/Rscript

# This is a script that takes as input a .csv file
# with the position, aminoacid and predicted disorder score (VSL2)
# and outputs a graph with the predicted disorder score along
# the length of the protein

# Load the required libraries
library("ggplot2")

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Check the number of arguments
if (length(args)!=1) {
	print("Please provide the correct number of arguments...")
	print("USAGE: Rscript plot_disorder_score.r [INPUT .csv]")
	q(save="no")
}

# Read the file and plot the graph
score <- read.csv(sprintf("%s", args[1]), header=TRUE)
colnames(score) <- c("position","AA","disorder")
head(score)

# Plot the line
output <- gsub(".csv",".pdf",args[1])
pdf(output)
ggplot(data=score, aes(x=position, y=disorder, group=1)) + geom_line() + xlab("Residue Number") + ylab("PONDR Score") + geom_line() + geom_hline(yintercept=0.5, color="red", linetype="dashed") + theme_bw() 
dev.off()
