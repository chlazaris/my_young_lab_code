#!/bin/Rscript

#This is a script to pick "random" samples

# Read the command line arguments
input <- commandArgs(trailingOnly=TRUE)

if(length(input)!=2){
	print("Please provide the correct number of arguments...")
	print("USAGE: Rscript random_number_generator.R [INPUT .tsv file] [Number of samples picked]")
	q(save="no")
}

# Read the file with the samples
data <- read.table(sprintf("%s", input[1]), header=TRUE, sep="\t")

# Get the total number of samples
total_no <- nrow(data)

# Get the number of picked samples
selected_no <- as.numeric(sprintf("%s", input[2]))

# Select the random samples (with no replacement)
selected_samples <- sample(1:total_no, selected_no, replace=FALSE)

# Generate the "number" column of the data
data$number <- c(1:total_no)

# Generate the "picked" (yes/no) column
data$picked <- ifelse(data$number %in% selected_samples, 'yes', 'no') 

# Get the output .tsv file
input_file <- sprintf("%s", input[1])
outfile <- gsub(".tsv", "_final.tsv", input_file)
write.table(data, outfile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
