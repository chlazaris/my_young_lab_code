#!/bin/bash

# Get the input .bam file
input=$1

# Check if the right number of arguments has been provided
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: create_bam_index [INPUT BAM FILE]"
	exit 1
fi

# Create the index
samtools index $input ${input}.bai
