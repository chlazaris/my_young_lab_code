#!/bin/bash

# Check for input
input=$1

if ($# != 1); then
	echo "Provide the right number of arguments..."
	echo "USAGE: convert_bam_to_sam.sh [INPUT BAM FILE]"
	exit
fi

# Convert .bam to .sam
outfile=${input%.bam}.sam
samtools view -h -o $outfile $input

# Compress the .sam file
gzip $outfile
