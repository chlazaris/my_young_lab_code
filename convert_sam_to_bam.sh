#!/bin/bash

# Provide input
input=$1

# Check for the number of arguments
if [ "$#" -ne 1 ]; then
	echo "Please provide input SAM file"
	echo "USAGE: ./convert_sam_to_bam INPUT_SAM"
fi

# Run samtools to convert SAM to BAM
samtools view -S -b $input > ${input%.sam}.bam 
