#!/bin/bash

if [ "$#" -ne 1 ]; then
	echo "Please provide the required arguments..."
	echo "USAGE: get_genebody.sh mm9, mm10, hg19, hg18..."
	exit
fi

# Get the input
genome=$1

# -N : no headers
    # -B : tab-delimted output
    # uniq to remove duplicate TSSs across tmultiple transcripts
    # grep -v "_" to remove unplaced contigs
    mysql --user genome \
          --host genome-mysql.cse.ucsc.edu \
          -N \
          -B \
          -D $genome \
          -e  "SELECT chrom, txStart, txEnd, \
                      X.geneSymbol, 1, strand \
               FROM knownGene as K, kgXref as X \
               WHERE txStart != txEnd \
               AND X.kgID = K.name" \
    | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' \
    | sort -k1,1 -k2,2n \
    | uniq \
    | grep -v "_" \
    | grep -vE 'random|^chrY|^chrM' \
    > ${genome}_knownGene.bed
