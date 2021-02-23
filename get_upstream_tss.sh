#!/bin/bash

# Specify input
if [ "$#" -ne 3 ]; then
	echo "Please provide the right number of arguments..."
	echo "USAGE: get_upstream_tss UPSTREAM_REGION_NT TSS_FILE CHROM_SIZES"
	exit
fi

# Specify inputs
upstream=$1
tss=$2
chrom_sizes=$3

# Specify the name of output
outfile=${tss%.bed}

bedtools slop \
	     -s \
             -l $upstream \
	     -r 0 \
             -i $tss \
             -g $chrom_sizes \
    > $outfile.upstream.${upstream}bp.bed
