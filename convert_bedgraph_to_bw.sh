#!/bin/bash

# Get the input
input=$1
output=${input%.bedGraph}.bw
chr_sizes=$2

# Check the number of inpus
if [ "$#" -ne 2 ]; then
	echo "Please provide 2 arguments"
	echo "USAGE: convert_bedgraph_to_bw BEDGRAPH CHROM_SIZES"
	exit
fi

# Convert to bigwig
ucsc_tools/bedGraphToBigWig $input $chr_sizes $output
