#!/bin/bash

# Get input .bam file
input=$1

# Create a bigwig directory
# mkdir bigwig

# Get the outfile
bedgraph_file=$(basename $input | sed 's/.bam/.BedGraph/')
#bedgraph_file2=$(basename $input | sed 's/.bam/.BedGraph2/')
bigwig_file=$(basename $input | sed 's/.bam/.bw/')

# Run bamcoverage with RPKM normalization
# bamCoverage -b $input --normalizeUsingRPKM -o bigwig/$outfile

# Convert to reads per million of mapped reads
TmpScale=$(bc <<< "scale=2;1000000/$(samtools view -F 4 -c $input)")
bedtools genomecov -ibam $input -bg -scale $TmpScale -g hg19.chrom.sizes | sortBed > bigwig/$bedgraph_file

# Change the bedgraph to work with chrom sizes
#cat bigwig/$bedgraph_file | grep -v MT | sed -i 's/^/chr/' > bigwig/$bedgraph_file2

# Convert .bedGraph to BigWig
bedGraphToBigWig bigwig/$bedgraph_file hg19.chrom.sizes bigwig/$bigwig_file

# Remove the useless files
#rm -rf bigwig/$bedgraph_file2
#rm -rf bigwig/$bedgraph_file
