#!/bin/bash -l

# Load the required modules
module load python/2.7

if [ $# -eq 3 ]; then
	continue
elif [ $# -eq 4 ]; then
	continue
else
	echo "Please provide the right number of arguments..."
	echo "USAGE: run_rose2 GFF-FILE SIGNAL-FILE INPUT-FILE [DISTANCE FROM TSS (to exclude)]"
        exit 1
fi

gff=$1
signal=$2 #.bam with signal (typically H3K27ac)
input=$3
tss_distance=${4:-0}

# Define output
output=$(echo $signal | cut -d'.' -f1)

# Run ROSE2 
rose2 -g HG19 -i $gff -r $signal -c $input $tss_distance -o ${output}_supers       
