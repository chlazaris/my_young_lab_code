#!/bin/bash

# Script to convert the input (3 column format) to SAF format 

# Get the input
input=$1
prefix=$2 # This is to be used to create a unique identifier

# Check for the number of arguments
if [ "$#" -ne 2 ]; then
	echo "Please provide right number of arguments..."
	echo "convert_3col_to_saf INPUT_FILE PREFIX"
	exit
fi

# Get the output name
out=$(basename $input)
outfile=${out%.bed}.saf

# Convert the input file to SAF format
echo "Creating the .SAF file..."
cat $input | awk -v prefix=$prefix '{print prefix"_"NR"\t"$1"\t"$2"\t"$3"\t""."}' > $outfile
echo "Done!"
