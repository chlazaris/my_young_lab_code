#!/bin/bash

# Get the peaks
peaks=$1

# Create a directory to save everything
# The directory should have the name of the peaks file
current_dir=$PWD
outdir=$(basename $peaks | cut -d'_' -f1)

if [ -d $outdir ]; then
        echo "The directory $outdir already exists..."
        exit
else
        mkdir $outdir
fi

# Go to the new directory
cd $outdir

# Create the links
ln -s /lab/solexa_young/lazaris/code


# Compute the matrix for each one of the .bw files
code/run_computeMatrix.sh $peaks

# Now generate the heatmap with the
# average profile above
code/run_plotHeatmap.sh descend_matrix.gz

# Get back to the previous directory
cd $current_dir
