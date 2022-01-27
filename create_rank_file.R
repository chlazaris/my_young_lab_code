#!/bin/Rscript

# This is a script that reads in an RNA-Seq result
# from the RNA-Seq pipeline and creates the corresponding
# rank (.rnk) file

# Get the arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1){print("USAGE: Rscript create_rank_file.R RNA-Seq-result.tsv"); quit(save="no")}

# Get the input
input=args[1]

# Read in the data
data <- read.table(sprintf("%s", input), header=TRUE, check.names=FALSE)
cols <- c(1,2,3)

# Select the columns you want
data2 <- data[,cols]

# Create the rank
data2$rank <- sign(data2$log2FoldChange) * (-log10(data2$pvalue))

# Save the ordered genes and the rank without header
data_rnk <- data2[order(data2$rank, decreasing=TRUE),]

# Get the output
out <- strsplit(input,"[.]")[[1]][1]
out1 <- paste(out,"rnk",sep=".")

# Now save the ranked file
cols <- c(1,4)
data_f <- data_rnk[,cols]
write.table(data_f,sprintf("%s",out1),col.names=F,row.names=F,quote=F,sep="\t")
