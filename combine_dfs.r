#!/bin/Rscript

# Read the files into a list
l <- lapply(Sys.glob("*final.featureCounts"), read.table, header=TRUE)

# Convert the list to dataframe
df <- do.call(rbind.data.frame, l)
df$length_in_kb <- (df$Length)/1000
df$RPK <- df$no_reads/df$length_in_kb
dim(df)

write.table(df,"all_data_for_enhancer_SE.tsv", col.names=T, row.names=F, sep="\t", quote=F)
