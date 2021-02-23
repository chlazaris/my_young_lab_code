#!/bin/Rscript

library("ggpubr")

# Get the input .csv file
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("Provide the right number of arguments..."); q(save="no")}
input <- args[1]

# Read in the data
data <- read.csv(sprintf("%s", input), header=T)

# Get the two columns as strings
col1 <- sprintf("%s", names(data)[1])
col2 <- sprintf("%s", names(data)[2])

# Get the name of output
output <- gsub("csv","pdf",input)
print(output)
print(data)

# Print the plot
pdf(output)
ggbarplot(data, x = col1, y = col2,
          fill = "blue",               # fill color
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45           # Rotate vertically x axis texts
          )
dev.off()
