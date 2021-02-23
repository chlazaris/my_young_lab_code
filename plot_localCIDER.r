#!/bin/Rscript

# Load the required libraries
library("ggplot2")

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)

# Check the number of arguments
if(length(args)!=1){
	print("Please provide the localCIDER .csv file")
	print("USAGE: Rscript plot_localCIDER.r [.csv]")
	q(save="no")
}

# Get the .csv file
cider_data <- read.csv(sprintf("%s", args[1]), header=T)
cider_data$AA <- c(1:nrow(cider_data))
cider_data$NCPR

# Create the output
output <- gsub(".csv","_plot.pdf",args[1])
pdf(output)
ggplot(data=cider_data, aes(x = AA, y = NCPR)) + 
       geom_area(data = subset(cider_data, cider_data$NCPR <= 0), fill = "red") +
       geom_area(data = subset(cider_data, cider_data$NCPR >= 0), fill = "blue") + ylim(-0.3,0.3) +
       scale_x_continuous(expand = c(0, 0)) + 
       theme_bw()
dev.off()
