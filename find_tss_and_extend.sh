#!/bin/bash

# This is a script that extracts the TSS
# and extends it on both sides
# by the specified number of basepairs

# Get the necessary inputs
input=$1
extend=$2
chr_sizes=$3

# Check for the required input
if [ $# -ne 3 ]; then
	echo "Please provide the required arguments..."
	echo "USAGE: find_tss_and_extend.sh [.bed with gene coordinates] [extension in bp] [chrom sizes]"
	exit 1
fi

# Specify the extension in bp
# to get it back in kilobases
extend_kb=$(echo $extend / 1000 | bc)

# Get name for the output file
gene_tss=${input%.coord.bed}.tss.bed

# Specify the TSS
# based on the strand
awk -F "\t" '{if ($6 == "+") {print  $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} \
	else {print  $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}  }' $input > $gene_tss

# Specify the name for the file with the extended TSS
extended_gene_tss=${gene_tss%.bed}.plus_minus_${extend_kb}kb.bed

# Extend the TSS taking into account the strand
slopBed -i $gene_tss -g $chr_sizes -l $extend -r $extend -s | sort | uniq | sortBed > $extended_gene_tss 
