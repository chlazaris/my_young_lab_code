#!/bin/bash

# Get the matrix file as input
matrix=$1

# Check for the number of args
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: plot_heatmap.sh MATRIX_FILE"
	exit 1
fi

# Get the name for the output heatmap file
outfile=${matrix%.gz}.pdf

# Plot a heatmap
plotHeatmap -m $matrix --dpi 300 --colorMap Blues --missingDataColor '#FFFFFF' -o $outfile
