#!/bin/bash

# This script is used to add chroms to bigwig files (from the Nextflow pipeline)
# and remove chrMT

input=$1

# Check for the number of arguments
if [ $# -ne 1 ]; then
	echo "Please provide the right number of arguments..."
	echo "USAGE: add_chr_to_bigwig.sh [INPUT bigWig file]"
	exit 1
fi

# Location of chrom sizes for hg19
hg19_path="/lab/solexa_young/lazaris/genomes/hg19"

# Convert bigwig to bedgraph
outfile=${input%.bigWig}.bdg
bigWigToBedGraph $input $outfile

# Replace chrom and remove MT chromosome
outfile2=${outfile%.bdg}_chr.bdg
cat $outfile | sed 's/^/chr/' | grep -v chrMT > $outfile2

# Convert the new BedGraph to BigWig again
bedGraphToBigWig $outfile2 $hg19_path/hg19.chrom.sizes ${outfile2%.bdg}.bigWig

