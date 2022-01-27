#!/bin/bash

# Run in on the cluster
eval "$(conda shell.bash hook)"

# activate the HOMER environment
conda activate homer

# Generate peaks starting with 'chr'
cat peaks/IR_peaks.narrowPeak | sed 's/^/chr/' > peaks/IR_peaks_chr.narrowPeak

# Annotate the peaks
annotatePeaks.pl peaks/IR_peaks_chr.narrowPeak hg19 > annotated_IR_peaks.tsv

# deactivate the HOMER environment
conda deactivate homer
