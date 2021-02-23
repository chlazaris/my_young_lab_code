#!/bin/Rscript

# Check for required packages
if (!require(dplyr)) install.packages('dplyr')

# Load required packages
library(dplyr)

# Specify the number of arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {print("Please provide the right arguments... USAGE: inner_join.r FILE1 FILE2 COMMON_FIELD OUTFILE"); q(save="no")}

# Read in the tables
x <- read.table(sprintf("%s", args[1]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
y <- read.table(sprintf("%s", args[2]), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Remove rows from x,y where something is missing
x[x==""] <- NA
y[y==""] <- NA
x <- na.omit(x)
y <- na.omit(y)

common_id <- sprintf("%s", args[3])
outfile <- sprintf("%s", args[4])

# Perform inner join
df <- inner_join(x, y, by = common_id)
write.table(df, outfile, row.names=F, col.names=T, quote=FALSE, sep="\t")
