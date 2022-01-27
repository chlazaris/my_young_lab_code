#!/bin/bash

# Specify the path to CRCmapper
crcmapper_path="/nfs/BaRC_Public/BaRC_code/Python/CRCmapper/"

# Get the required inputs
enhancers=$1
input_bam=$2
fasta=$3
subpeaks=$4

# Run CRCmapper
python2 $crcmapper_path/CRCmapper.py -e input/$enhancers -b input/$input_bam -g HG19 -f input/$fasta -s input/$subpeaks -n MM1S_CRCmapper_CRC -o MM1S_CRCmapper_CRC  
