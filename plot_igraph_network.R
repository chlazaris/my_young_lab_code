#!/bin/Rscript

library(igraph)

# Specify the input
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]

# Read in the file with the edges
data <- read.table(sprintf("%s", input), header=T, sep="\t", stringsAsFactors=FALSE)

# Convert to edge matrix
m <- as.matrix(data)

# Create the output file name
outfile <- gsub("CRC_edges.tsv","network.pdf",input)

# Create and graph a directed network
network <- graph_from_edgelist(m, directed=TRUE)

pdf(outfile, useDingbats=FALSE)
plot(network, layout=layout.circle, vertex.color=rgb(0.255,0.412,0.882), vertex.label.family="Helvetica", vertex.label.font=7, 
     vertex.label.color="white", vertex.size=32, edge.width=0.5, edge.arrow.size=0.5, edge.arrow.width=0.4, edge.lty=c("dashed"))
dev.off()
