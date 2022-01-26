#!/bin/bash

# This is a script that takes as input 
# the *AllEnhancers.table.txt from ROSE2
# and outputs .bed files for
# a. all enhancers (AllEnhancers.bed)
# b. typical enhancers (*TypicalEnhancers.bed)
# c. super-enhancers (*SuperEnhancers.bed)

# Check for input
if [ $# -ne 1 ]; then
	echo "Please, provide the correct number of arguments..."
	echo "USAGE: extract_enhancers_from_table.sh [INPUT; *AllEnhancers.table.txt]"
	exit 1
fi

input=$1

all_enhancer=${input%.table.txt}.bed

# Get all enhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
	awk '{print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
	sort | uniq | sortBed > $all_enhancer

# Get the number of all enhancers
echo "All enhancers: $(wc -l $all_enhancer | cut -d' ' -f1)"

typical_enhancer=$(echo $all_enhancer | sed 's/All/Typical/')

# Get the typical enhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
	awk '{ if ($10==0) print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
        sort | uniq | sortBed > $typical_enhancer

# Get the number of typical enhancers
echo "Typical enhancers: $(wc -l $typical_enhancer | cut -d' ' -f1)"

super_enhancer=$(echo $all_enhancer | sed 's/All/Super/')

# Get the superenhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
        awk '{ if ($10==1) print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
        sort | uniq | sortBed > $super_enhancer

# Get the number of typical enhancers
echo "Super-enhancers: $(wc -l $super_enhancer | cut -d' ' -f1)"
