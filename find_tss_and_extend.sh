#!/bin/bash

#This is a script that extract the TSSs 

input=$1
extend=$2

extend_kb=$(echo $extend / 1000 | bc)

awk -F "\t" '{if ($6 == "+") {print  $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} else {print  $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}  }' $input > ${input%.bed}.tss.bed

# Get name for the output file
gene_tss=${input%.bed}.tss.bed

extended_gene_tss=${gene_tss%.bed}.plus_minus_${extend_kb}kb.bed
# Get a certain window around the TSS
slopBed -i $gene_tss -g hg19.chrom.sizes -l $extend -r $extend -s | sort | uniq | sortBed > $extended_gene_tss 
