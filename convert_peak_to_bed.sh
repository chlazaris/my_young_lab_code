#!/bin/bash

# Get the input
input=$1

# Get the relevant columns
cut -f1-6 $input | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3,$4,"0",$6}' | sortBed > ${input%.narrowPeak}.bed
