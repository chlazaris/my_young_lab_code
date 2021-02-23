#!/bin/bash

# Get the inputs
features=$1
bam_file=$2
output=$3

# Calculate the per-base coverage over certain features
bedtools coverage -a $features -b $bam_file -d > $output  
