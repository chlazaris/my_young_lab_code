#!/bin/bash

# Give the .bam file as input
input=$1
output=${input%_sorted.bam}

# Check the number of arguments
if [ "$#" -ne 1 ]; then
	echo "Please provide correct number of arguments..."
	echo "USAGE: ./map_to_strand SORTED_BAM"
	exit
fi

# Get the reads per strand
bedtools genomecov -ibam $input -bg -strand + > ${output}_plus.bedGraph
bedtools genomecov -ibam $input -bg -strand - > ${output}_minus.bedGraph
