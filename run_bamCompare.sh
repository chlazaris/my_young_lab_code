#!/bin/bash

# Check for the number of arguments
if [ $# -ne 2 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: bamCompare [SIGNAL .bam] [INPUT .bam]"
	exit
fi

# Get the signal and the input
signal=$1
input=$2

# Subtract the input from the signal file
# and normalize based on BPM (bins per million)
outfile=$(basename $signal | sed 's/.bam//')_BPM.bw
#echo $outfile

# Subtract the input and normalize using BPM
bamCompare -b1 $signal -b2 $input --normalizeUsing BPM --operation subtract --scaleFactorsMethod None --effectiveGenomeSize 2864785220 -o $outfile    
