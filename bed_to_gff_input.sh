#$/bin/bash

# Get the input .bed file
input=$1
category=$2

# Convert the .bed file for input to script for .gff conversion
cat $input | awk -v category="$category" 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4"_"category"_"NR,".",$6}' | sortBed > ${input%.bed}_f.bed
