#!/bin/Rscript

# Load the required libraries
library(ENCODExplorer)

# Get the snapshot 
encode_df <- get_encode_df()

# Get the results you care about
query_results <- queryEncode(organism = "Homo sapiens", target=c("H3K27ac","H3K4me1","H3K4me3","Control"),
                      biosample_name = "natural killer", file_format = "fastq",
                      fixed = TRUE, fuzzy=TRUE)
# Write the results to a file
result_df <- data.frame(query_results)
write.table(result_df, "NK_datasets_from_ENCODE.txt", row.names=F, col.names=T, sep="\t", quote=F)
