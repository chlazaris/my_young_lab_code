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
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: run_gimmemotifs.sh [INPUT-FILE (.bed with peak summits)]"
	exit 1
fi

# Get the output directory
outdir=${input%.bed}.motifs

# Run gimmemotifs
gimme motifs $input $outdir -p JASPAR2020_vertebrates -g hg19 --known

# Deactivate the conda environment
conda deactivate
