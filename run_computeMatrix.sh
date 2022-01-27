#!/bin/bash

# Get the number of arguments
#if [ "$#" -ne 2 || "$#" -ne 3 || "$#" -ne 4 ]; then
#	echo "Please provide the correct number of arguments..."
#	echo "USAGE: run_computeMatrix.sh [.bed peaks] [file with .bw] [extend; default: 2kb] [order; default: descend]"
#	exit
#fi	

# Get the input (sort type) or set it to default
peaks=$1
file=$2
extend=${3:-2000}
order=${4:-descend}

# Create a directory to save everything
# The directory should have the name of the peaks file
current_dir=$PWD
outdir=$(basename $peaks | cut -d'_' -f1)

if [ -d $outdir ]; then
	echo "The directory $outdir already exists..."
	exit
else
	mkdir $outdir
fi

#Generate the file with the relative paths to .bw
cat $file | sed 's/^/bigwig\//' | sed 's/$/.bw/' > $outdir/bigwig.txt

# Go to the new directory to generate the matrix
cd $outdir

# Create the necessary links
ln -s ../peaks
ln -s ../../bigwig
bw_files=$(cat bigwig.txt | tr '\n' '	')

# Compute the matrix
computeMatrix reference-point --referencePoint center -S $bw_files \
	-R $peaks -bs 50 --sortRegions $order --sortUsing mean \
	-a $extend -b $extend -p 8 --skipZeros -o ${order}_matrix.gz \
        --outFileNameMatrix ${order}IndividualValues.tsv --outFileSortedRegions ${order}SortedRegions.bed
