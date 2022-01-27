#!/bin/bash

# This is a script that runs gimmemotifs on the cluster
# using hg19 as the reference genome and JASPAR2020 as
# the motif database (for known motifs)

# Activate the environment
eval "$(conda shell.bash hook)"
conda activate gimmemotifs0_16_1

# Set cache properly
NEW_CACHE=./cache
mkdir -p $NEW_CACHE
if [ -z $XDG_CACHE_HOME ]; then
    XDG_CACHE_HOME=$HOME/.cache
fi
cp -r $XDG_CACHE_HOME/gimmemotifs $NEW_CACHE/
export XDG_CACHE_HOME=$NEW_CACHE
echo 'Using $XDG_CACHE_HOME for cache'

# Check the number of arguments
if [ $# -ne 2 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: run_gimmemotifs.sh [INPUT-FILE (peak file; narrowPeak, bed)] [MOTIF status (known, denovo)]"
	exit 1
fi

# Get the input
input=$1
motif_status=$2

# Get the output directory
outdir=$(basename $input | sed 's/bed/motifs/')

# Run gimmemotifs
if [ "$motif_status" == "known" ]; then
        gimme motifs $input $outdir -p JASPAR2020_vertebrates -g hg19 --size 200 --known
elif [ "$motif_status" == 'denovo' ]; then
        gimme motifs $input $outdir -g hg19 --size 200 --denovo
else
	echo "Invalid motif status: Please provide motif status known or denovo"
fi

# Deactivate the conda environment
conda deactivate
