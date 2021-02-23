#!/bin/Rscript

# This is a script to graph disk usage (top 10 users)
# Gives the output in TB

# Load the required libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library("ggplot2")
library("dplyr")

# Specify the command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {print("Provide correct number of args. USAGE: graph_disk_space CSV_file"); q(save="no")}

# Get the input name
input <- args[1]

# Read in the .csv file
data <- read.csv(sprintf("%s", args[1]), header=T, stringsAsFactors=F, check.names=F)
data2 <- data %>% select(UID, Bytes) %>% data.frame()
data3 <- data2 %>% mutate(storage_in_TB=Bytes/10^12) %>% data.frame()
ranked <- data3[order(-data3$storage_in_TB),]
data_to_plot <- ranked[1:10,]

# Get the output name
out <- gsub(".csv", ".pdf", input)

pdf(out)
ggplot(data_to_plot, aes(x=reorder(UID, storage_in_TB), y=storage_in_TB)) + geom_bar(stat="identity", fill="#0000FF") + coord_flip() + xlab("User ID") + ylab("Storage usage (in TB)") + theme_bw()
dev.off()
