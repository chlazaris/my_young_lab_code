#!/bin/bash

# This is a script that takes as input the total peaks, enhancer peaks and promoter
# peaks and calculates their numbers and percentages

if [ $# -ne 3 ]; then
	echo "Please provide the correct number of arguments..."
	echo "USAGE: peak_stats.sh [TOTAL PEAKS] [PROMOTER] [ENHANCER]"
	exit 1
fi

# Get the numbers
total=$(cat $1 | wc -l)
promoter=$(cat $2 | wc -l)
enhancer=$(cat $3 | wc -l)
echo $total
echo $promoter
echo $enhancer

# Get fractions
prom=$(echo "scale=2; $promoter/$total" | bc -l)
enh=$(echo "scale=2; $enhancer/$total" | bc -l)
other=$(echo "scale=2; 1-$prom-$enh" | bc -l)

# Get the results in a peak stats file
if [ ! -f peak_stats.txt ]; then
	echo "Total_peaks: $total" > peak_stats.txt
        echo "Promoter_peaks: $promoter" >> peak_stats.txt
	echo "Enhancer_peaks: $enhancer" >> peak_stats.txt
	echo "Promoter_fraction: $prom" >> peak_stats.txt
	echo "Enhancer_fraction: $enh" >> peak_stats.txt
	echo "Other_fraction: $other" >> peak_stats.txt
else
	echo "File with peak stats already exists..."
	exit
fi
