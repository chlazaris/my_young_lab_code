#!/bin/bash

# This is a script used to clean up the output from ENCODEExplorer and keep only the relevant
# columns

# Check for the right number of arguments
if [ $# -ne 1 ]; then
	echo "Please provide the ENCODE explorer output as input"
	echo "USAGE: cleanup_encode.sh FILE.txt"
	exit 1
fi

# Get the input file
input=$1

# Clean up the file
cat $input | cut -f8,11,22,23,28,31,37,38,48,49,51,54 > ${input%.txt}_clean.txt
