#!/bin/bash

# Get the input (fasta file for transcripts)
input=$1

# Check the input
if [ "$#" -ne 1 ]; then
	echo "Please provide the right number of arguments..."
	echo "./build_salmon_index FASTA"
	exit
fi

# Generate the index
salmon index -t $input -i hg19_transcripts_index
