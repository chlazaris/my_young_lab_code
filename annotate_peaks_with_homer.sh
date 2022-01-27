#!/bin/bash

# Activate the environment
eval "$(conda shell.bash hook)"
conda activate homer

# Get the peaks as input (bed format including narrowPeak)
input=$1
input_name=$(echo $input | cut -d'.' -f1)

# Perform the annotation
annotatePeaks.pl $input hg19 > ${input}.annotated

# Deactivate the environment
conda deactivate
