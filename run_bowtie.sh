#!/bin/bash

# Get the input (fastq file)
input=$1

# Check the number of arguments
if [ "$#" -ne 1 ]; then
	echo "Please provide right number of arguments..."
	echo "USAGE: ./run_bowtie FASTQ_FILE"
	exit
fi

# Run Bowtie
bowtie --index /lab/solexa_young/lazaris/genomes/Mus_musculus/references/UCSC/mm9/Sequence/BowtieIndex/genome -e 70 -k 1 -m 10 -n 2 --best -q $input -S ${input%.fastq}.sam
