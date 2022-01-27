#!/bin/bash

# Check for the number of arguments
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: bamCoverage [SIGNAL .bam]"
	exit
fi

# Get the signal
signal=$1

# Set the path to blacklist path
blacklist_path="/lab/solexa_young/lazaris/genomes/hg19/blacklist/"

# Subtract the input from the signal file
# and normalize based on BPM (bins per million)
outfile=$(basename $signal | sed 's/.bam//').bw
#echo $outfile

# Subtract the input and normalize using BPM
bamCoverage -b $signal \
	--binSize 50 \
	--normalizeUsing CPM \
	--extendReads 150 \
	--ignoreDuplicates \
        --ignoreForNormalization chrX chrY chrM \
	--effectiveGenomeSize 2620345972 \
	--blackListFileName $blacklist_path/hg19_blacklist_no.bed \
	--numberOfProcessors 6 \
	-o bigwig/$outfile
