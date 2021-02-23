#!/bin/Rscript

# Plot polar plot

# Install the required packages
if (!require("fmsb")) install.packages("fmsb")

# Load the required libraries
library(fmsb)

# Read in the arguments
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]

# Read in the data
data <- read.csv(sprintf("%s", input), header=TRUE, stringsAsFactors=F)
data$logFDR <- -log(data$FDR)
order.data <- order(data$logFDR, decreasing=T)
ordered_data <- data[order.data,]
ordered_data <- ordered_data[1:10,]

# Convert from long to wide
data2 <- as.data.frame(matrix(ordered_data$logFDR, nrow=1, ncol=nrow(ordered_data)))
colnames(data2) <- ordered_data$term_name

# Find the number of categories 
no_categories <- length(ordered_data$term_name)

# Find the max and min of logFDR
max_logFDR <- max(data$logFDR)
min_logFDR <- min(data$logFDR)

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min 
data3 <- rbind(rep(max_logFDR, no_categories) , rep(min_logFDR, no_categories) , data2)

# Check your data, it has to look like this!
head(data3)

# Custom the radarChart !
pdf("test_polar.pdf", width=15, height=10)
radarchart(data3, axistype=1 ,
    #custom polygon
    pcol=rgb(0.447,0.737,0.831) , pfcol=rgb(0.678,0.847,0.902) , plwd=2 , pdensity=0.5,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
    #custom labels
    vlcex=0.8
    )
dev.off()
