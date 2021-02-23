#!/bin/bash

if [ "$#" -ne 1 ]; then
       echo "Incorrect number of arguments"
       echo "USAGE: ./sort_index_bam.sh INPUT_BAM"       
fi

# Provide the .bam file as input
input=$1
output=${input%.bam}_sorted.bam

# Sort and index the bam file
samtools sort $input > $output
samtools index $output
