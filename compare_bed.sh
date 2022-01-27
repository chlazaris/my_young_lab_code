#!/bin/bash

# This is a script that compares
# two .bed files (A and B) and 
# returns:
# a. The items in A, not in B
# b. The items in A that are also in B

file1=$1
file2=$2
option=$3

if [ $# -ne 3 ]; then
	echo "USAGE: compare_bed.sh [.bed FILE 1] [.bed FILE 2] [option: common or diff]"
	exit
fi

# Check the type of comparison and perform it
if [ $option = "common" ]; then
	outfile1=$(basename $file1 | cut -d'.' -f1)"_in_"$(basename $file2 | cut -d'.' -f1)".bed"
	bedtools intersect -wa -a $file1 -b $file2 | sort | uniq | sortBed > $outfile1
elif [ $option = "diff" ]; then
	outfile2=$(basename $file1 | cut -d'.' -f1)"_not_in_"$(basename $file2 | cut -d'.' -f1)".bed"
	bedtools intersect -wa -a $file1 -b $file2 -v | sort | uniq | sortBed > $outfile2
else
	echo "Please, provide a valid option..."
	exit 1
fi	
