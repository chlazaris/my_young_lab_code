#!/bin/bash

# Get the input files
saf_file=$1
bam_file=$2

# Check if the input is correct
if [ "$#" -ne 2 ]; then
	echo "Please provide the right number of arguments..."
	echo "USAGE: run_featureCounts INPUT_SAF INPUT_BAM"
	exit
fi

# Run featureCounts
echo "Running..."
featureCounts --primary -F 'SAF' -a $saf_file -s 0 -o ${saf_file%.saf}.featureCounts $bam_file
echo "Done!"
