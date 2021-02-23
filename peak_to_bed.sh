#!/bin/bash

if [ $# -ne 1 ]; then
	echo "Please provide the right number of arguments..."
	echo "peak_to_bed INPUT_FILE (e.g. narrowPeak)"
	exit
fi

# Specify the input file
input=$1

# Convert to .bed
cut -f1-6 $input | sed 's/^/chr/' | sortBed > ${input%.narrowPeak}.bed
