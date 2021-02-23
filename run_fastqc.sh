#!/bin/bash

# Specify input
input=$1
outdir=$2

# Test number of arguments
if [ "$#" -ne 2 ]; then
    echo "Illegal number of input files"
    echo "USAGE: run_fastqc OUTDIR INPUT_FASTQ"
    exit
fi

#Create directory for the results
mkdir $outdir

# Run fastqc and send results to the output directory
fastqc -o $outdir $input
