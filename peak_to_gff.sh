#!/bin/bash

# Get the input narrowPeak or broadPeak file
input=$1

# Check if the number of arguments is correct
if [ $# -ne 1 ]; then
	echo "Please provide the right number of arguments..."
	echo "USAGE: peak_to_gff [INPUT narrowPeak/broadPeak file]"
	exit 1
fi	

# Convert it to gff
awk -F\\t '{print $1 "\t" $4 "\t\t" $2 "\t" $3 "\t\t.\t\t" $4}' $input > ${input}.gff.tmp
sed 's/^/chr/' ${input}.gff.tmp > ${input}.gff
rm -rf ${input}.gff.tmp
