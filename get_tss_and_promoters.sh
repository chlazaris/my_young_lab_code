#!/bin/bash

# This script extracts the gene coordinates using
# Refseq genes as input (the .gtf file from UCSC; hg19.refGene.gtf)

# Check for the number of arguments
if [ $# -ne 2 ]; then
	echo "Please provide the correct number of arguments..."
        echo "USAGE get_tss_and_promoters.sh [extension] [chromosome sizes]"
        exit 1
fi

# Get the input files
extension=$1
chr_sizes=$2

# Specify the name of the .gtf file
gtf_file="hg19.refGene"

# Download the .gtf if it does not exist
if [ ! -f ${gtf_file}.gtf.gz ]; then
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz" -O ${gtf_file}.gtf.gz
fi

# Check if the decompressed file exists
if [ ! -f ${gtf_file}.gtf ]; then
	zcat ${gtf_file}.gtf.gz > ${gtf_file}.gtf
fi

# Extract the transcript coordinates
# Remove chromosomes M and Y (as Y not everywhere)
if [ ! -f ${gtf_file}.bed ]; then
	cat ${gtf_file}.gtf | gtf2bed | grep -w transcript | grep -vE 'Un_|hap|random|chrM|chrY' | cut -f1-6 | sort | uniq | sortBed > ${gtf_file}.bed
fi

# Extract the TSS
if [ ! -f ${gtf_file}.tss.bed ]; then
	cat ${gtf_file}.bed | awk '{if ($6 == "+") {print  $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} else {print  $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}}' | sort | uniq | sortBed > ${gtf_file}.tss.tmp
	bedtools merge -s -c 4 -o distinct -i ${gtf_file}.tss.tmp \
		| awk '{print $1"\t"$2"\t"$3"\t"$5"\t"".""\t"$4}' | sort | uniq | sortBed > ${gtf_file}.tss.bed
	rm -rf ${gtf_file}.tss.tmp	
fi

# Extend the TSS
ext_in_kb=$(echo $extension/1000 | bc)
if [ ! -f ${gtf_file}.tss.plus.minus.${extension}kb.bed ]; then
	slopBed -i ${gtf_file}.tss.bed -s -g $chr_sizes -l $extension -r $extension | sort | uniq | sortBed > ${gtf_file}.tss.plus.minus.${ext_in_kb}kb.bed	
fi

# Merge transcripts or genes that
# overlap
# to get a TSS per gene
#cat $out1 | mergeBed -s -d 0 \
#	-c 4 -o distinct \
#	| awk '{print $1"\t"$2"\t"$3"\t"$5"\t"".""\t"$4}' \
#	| sort | uniq | sortBed > ${out1%.bed}.coord.bed
