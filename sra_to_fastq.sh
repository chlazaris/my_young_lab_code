#!/bin/bash

input=$1

# Convert SRR to fastq (no split paired reads)
#fasterq-dump $input
#fasterq-dump $input --split-files
fasterq-dump $input

# Convert the fastq file to the .gz version
gzip ${input}.fastq

# Rename the file
#mv ${input}.fastq.gz ${input}_R1.fastq.gz
