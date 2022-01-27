#!/bin/bash

# Get the number of arguments
#if [ "$#" -ne 2 || "$#" -ne 3 || "$#" -ne 4 ]; then
#	echo "Please provide the correct number of arguments..."
#	echo "USAGE: run_computeMatrix.sh [.bed peaks] [extend; default: 2kb] [order; default: descend]"
#	exit
#fi

# Create the necessary links
ln -s ../peaks
ln -s ../../bigwig
ln -s ../bigwig.txt

# Get the input (sort type) or set it to default
peaks=$1
extend=${3:-2000}
order=${4:-descend}

#Generate the file with the relative paths to .bw
bw_files=$(cat bigwig.txt | sed 's/^/bigwig\//' | sed 's/$/.bw/' | tr '\n' ' ')

# Compute the matrix
computeMatrix reference-point --referencePoint center -S $bw_files \
	-R $peaks -bs 50 --sortRegions $order --sortUsing mean \
	-a $extend -b $extend -p 8 --skipZeros -o ${order}_matrix.gz \
        --outFileNameMatrix ${order}IndividualValues.tsv --outFileSortedRegions ${order}SortedRegions.bed
