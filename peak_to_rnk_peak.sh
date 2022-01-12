#!/bin/bash

# Convert a narrowPeak or broadPeak file to a  where the peaks 
# are ranked based on the ChIP signal

input=$1
input_suffix=$(echo $input | cut -d'.' -f2)

if [ $# -ne 1 ]; then
	echo "Please provide the right number of arguments..."
	echo "USAGE: peak_to_rnk_peak.sh [narrowPeak or broadPeak file]"
	exit
fi

# Get the ranked .bed file (exclude also the M and Y chromosomes)
outfile=$(echo $input | cut -d'.' -f1).rnk.$input_suffix
cat $input | sort -k7,7nr | sed 's/^/chr/' | grep -Ev 'chrY|chrM' > $outfile
