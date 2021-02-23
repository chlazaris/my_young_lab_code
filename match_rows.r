#!/bin/Rscript
library(plyr)

args <- commandArgs(trailingOnly=TRUE)

# Get the dataframes
df1 <- read.table(sprintf("%s", args[1]), header=T, sep="\t")
df2 <- read.table(sprintf("%s", args[2]), header=T, sep="\t")

data <- match_df(df1, df2, on='Geneid')
dim(data)
