#!/bin/bash

# Input the bigWig file
input=$1

# Convert .bw to .bdg
bigWigToBedGraph $input ${input%.bw}.bdg

# Compress the file
bgzip ${input%.bw}.bdg

# Create tabix index
tabix -p bed ${input%.bw}.bdg.gz
