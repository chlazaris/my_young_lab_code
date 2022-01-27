#!/bin/bash

# Get the number of arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: run_computeMatrix.sh [.bigwig file(s)] [.bed peak file(s)] [sort_sample; default:all] [extend; default:2000]"
        echo "sort_sample: sorting based on this sample or all (default)"
	exit 1
fi

# Get the input
bigwig=$1
peaks=$2
sort_sample=$3
extend=${4:-2000}

if [ "$sort_sample" = "" ]; then 
	# Compute the matrix
	computeMatrix reference-point --referencePoint center -S $bigwig \
		-R $peaks -bs 10 --sortRegions descend --sortUsing mean \
		-a $extend -b $extend -p 8 --skipZeros \
		-o matrix.gz \
        	--outFileNameMatrix matrixIndividualValues.tsv \
		--outFileSortedRegions matrixSortedRegions.bed
else
	# Compute the matrix that is know sorted based on the first sample
	computeMatrix reference-point --referencePoint center -S $bigwig \
                -R $peaks -bs 10 --sortRegions descend --sortUsing mean \
                -a $extend -b $extend -p 8 --skipZeros --sortUsingSamples $sort_sample\
                -o sample${sort_sample}SortedMatrix.gz \
                --outFileNameMatrix sample${sort_sample}SortedIndividualValues.tsv \
                --outFileSortedRegions sample${sort_sample}SortedRegions.bed
fi
