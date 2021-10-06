#!/bin/bash

# Check the number of arguments
if [ $# -ne 1 ]; then
	echo "Please provide the correct number of arguments..."
       	echo "USAGE: motif_logo.sh [JASPAR ID FOR PROTEIN OF INTEREST]"
	exit 1
fi

# Activate the gimmemotifs environment
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

input=$1 # Enter the protein of interest, so that you plot the motif

# Plot the motif
gimme logo -p JASPAR2020_vertebrates -i $input -k information

# Deactivate the conda environment
conda deactivate
