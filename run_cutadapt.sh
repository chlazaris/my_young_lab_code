#!/bin/bash

# Specify the input
input=$1

# Check for the number of arguments
if [ "$#" -ne 1 ]; then
	echo "Please provide one .fastq file as input"
	echo "USAGE: cutadapt INPUT_fastq"
	exit
fi

cutadapt -m 24 -q 33 -a TCGTATGCCGTCTTCTGCTTG -o ${input%.fastq}_trimmed_hq.fastq $input
