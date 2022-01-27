#!/bin/bash

# This script takes a compressed matrix (matrix.gz)
# from DeepTools as input and also color (equal to the number of samples)
# and outputs the corresponding profile

# Input
matrix=$1
color=$2

if [ $# -ne 2 ]; then
	echo "Please provide the correct number of arguments..."       
	echo "USAGE: run_plotProfile.sh [INPUT: matrix.gz] [color; e.g., red]"
	exit
fi

# Plot the profile
plotProfile -m $matrix \
        --outFileSortedRegions profileSortedRegions.bed \
        --outFileNameData profileIndividualValues.tsv \
        --perGroup \
        --colors $color \
        --plotHeight 9 \
        --plotWidth 10 \
        -o profile.pdf \
        --dpi 300
