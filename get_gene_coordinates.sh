#!/bin/bash

# This script extracts the gene coordinates using
# Refseq genes as input (the .gtf file from UCSC; hg19.refGene.gtf)

# Check for the number of arguments
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
        echo "USAGE get_gene_coordinates.sh [INPUT e.g., refGene file from UCSC; hg19.refGene.gtf]"
        exit 1
fi

gtf_file=$1

# Extract the gene coordinates
out1=${gtf_file%.gtf}.bed
cat $gtf_file | grep -w transcript | \
	grep -vE "hap|Un_gl|chrY|chrM|random" | \
	awk '{print $1"\t"$4"\t"$5"\t"$14"\t"$6"\t"$7}' | \
	sed 's/"//g' | sed 's/;//' | \
	sort | uniq | sortBed > $out1

# Merge transcripts or genes that
# overlap
# to get a TSS per gene
# Also, exclude everything that does not
# belong to chromosomes 1..22, chrX
cat $out1 | mergeBed -s -d 0 \
	-c 4 -o collapse \
	| awk '{print $1"\t"$2"\t"$3"\t"$5"\t"".""\t"$4}' \
	| sort | uniq | sortBed > ${out1%.bed}.coord.bed
