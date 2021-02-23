#!/bin/bash

# Check the number of inputs
if [ $# -ne 1 ]; then
	echo "Please provide the input .fasta file..."
	echo "USAGE: run_barcode_plot.sh [FASTA]"
	exit 1
fi

# Get the input .fasta file
fasta=$1

# Create the output file
barcode=${fasta%.fa}_barcode_plot.pdf

Rscript code/barcode_plot.r $fasta $barcode
