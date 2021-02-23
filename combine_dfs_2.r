#!/bin/Rscript

# Read the files into a list
l <- lapply(Sys.glob("*charge*"), read.table, header=TRUE)

# Convert the list to dataframe
df <- do.call(rbind.data.frame, l)

# Write into output
write.table(df,"all_data_for_enhancer_SE.tsv", col.names=T, row.names=F, sep="\t", quote=F)
