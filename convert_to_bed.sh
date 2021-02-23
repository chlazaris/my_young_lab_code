#!/bin/bash

# Get the input
input=$1

# Convert the input to .bed (6 column)
cat $input | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"0"\t""."}' | sortBed > ${input}.bed
