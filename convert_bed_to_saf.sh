#!/bin/bash

# Convert .bed to .saf
input=$1

# Check the number of arguments
if [ "$#" -ne 1 ]; then 
	echo "Please provide the right number of arguments..."
	echo "USAGE: convert_bed_to_saf INPUT_BED"
	exit
fi

# Make the conversion 
echo "Converting..."
cat $input | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$6}' > ${input%.bed}.saf
echo "Done!"
