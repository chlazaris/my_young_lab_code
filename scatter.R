#!/bin/Rscript

#library(devtools)
#devtools::install_github('cttobin/ggthemr')

# Load the required libraries
library(ggthemr)
library(ggplot2)
library(ggrepel)

# Select theme
ggthemr('dust')

args <- commandArgs(trailingOnly=TRUE)

# Get the input file
input <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t")
output <- gsub("txt", "", args[1])
output
outname <- paste0(output,"pdf")

pdf(outname, useDingbats=FALSE)
p <- ggplot(input, aes(In_Degree, Out_Degree, label=Tf)) +
geom_point(color="black") + geom_text_repel() + xlab("in-degree") +
ylab("out-degree")
plot(p)
dev.off()
#Reset the theme
ggthemr_reset()
