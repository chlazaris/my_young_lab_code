#!/bin/bash

# This is a script that takes as input 
# the *AllEnhancers.table.txt from ROSE2
# and outputs .bed files for
# a. all enhancers (AllEnhancers.bed)
# b. typical enhancers (*TypicalEnhancers.bed)
# c. super-enhancers (*SuperEnhancers.bed)
# It also outputs a stats.txt file with
# the number of promoters, enhancers, SEs

# Check for input
if [ $# -ne 2 ]; then
	echo "Please, provide the correct number of arguments..."
	echo "USAGE: extract_enhancers_from_table.sh [INPUT; *AllEnhancers.table.txt] [INPUT; promoters.bed]"
	exit 1
fi

# Get the required inputs
input=$1
promoters=$2

# Get the number of promoters
echo "Promoters: $(wc -l $promoters | cut -d' ' -f1)" >> stats.txt

all_enhancer=${input%.table.txt}.bed

# Get all enhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
	awk '{print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
	sort | uniq | sortBed > $all_enhancer.temp

# Now exclude the promoters to get the actual enhancers
bedtools intersect -v -wa -a $all_enhancer.temp -b promoters.bed > $all_enhancer
rm -rf $all_enhancer.temp

# Get the number of all enhancers
echo "All enhancers: $(wc -l $all_enhancer | cut -d' ' -f1)" >> stats.txt

typical_enhancer=$(echo $all_enhancer | sed 's/All/Typical/')

# Get the typical enhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
	awk '{ if ($10==0) print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
        sort | uniq | sortBed > $typical_enhancer.temp

# Now exclude the promoters to get the actual typical enhancers
bedtools intersect -v -wa -a $typical_enhancer.temp -b promoters.bed > $typical_enhancer
rm -rf $typical_enhancer.temp

# Get the number of typical enhancers
echo "Typical enhancers: $(wc -l $typical_enhancer | cut -d' ' -f1)" >> stats.txt

super_enhancer=$(echo $all_enhancer | sed 's/All/Super/')

# Get the superenhancers in a bed file
cat $input | grep -v '#' | sed 1d | \
        awk '{ if ($10==1) print $2"\t"$3"\t"$4"\t"$1"\t""0""\t""."}' | \
        sort | uniq | sortBed > $super_enhancer.temp

# Now exclude the promoters to get the actual super-enhancers
bedtools intersect -v -wa -a $super_enhancer.temp -b promoters.bed > $super_enhancer
rm -rf $super_enhancer.temp

# Get the number of super-enhancers
echo "Super-enhancers: $(wc -l $super_enhancer | cut -d' ' -f1)" >> stats.txt
