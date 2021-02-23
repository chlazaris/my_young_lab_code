#!/bin/bash

code_dir="/lab/solexa_young/lazaris/code/"
file1=$(basename *_AllEnhancers.table.txt*)
file2=$(basename *SuperEnhancers_ENHANCER_TO_TOP_GENE.txt)

# Remove the header so that you can join
file3=${file1%.table.txt}_final_table.txt
cat $file1 | grep -v "#" > $file3

# Perform the left join, to get everything in one file
Rscript $code_dir/left_join.r $file3  $file2 REGION_ID test.txt

# Create the input for the swoosh plot
cat test.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$19,$20,$21,$24,$25}' FS='\t' OFS='\t' > swoosh_plot_input.tsv

# Plot the super-enhancer plot
Rscript $code_dir/plot_superenhancer.r swoosh_plot_input.tsv

# Remove the test file as it is not needed anymore
rm -rf test.txt
