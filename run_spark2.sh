#!/bin/bash

# Get the current directory
mydir=$pwd

# Get into the other dir
#cd input

# Get the required inputs
coords=$1
gene=$2

#Code used to generate this plot:
python /lab/solexa_young/lazaris/tools/SparK/SparK.py \
-pr $coords \
-cf bedgraph/Donor1000_H3K27ac_R1.bdg bedgraph/Donor1022_H3K27ac_R1.bdg \
-gtf gencode.v19.annotation.gtf \
-gl Donor1000_H3K27ac Donor1022_H3K27ac \
-gs yes \
-sm 10 \
-bed bed/Donor1000_SuperEnhancers.bed bed/Donor1022_SuperEnhancers.bed \
-f 92C5DE 92C5DE \
-bedcol 000000 000000 \
-bedlab Donor1000_SE Donor1022_SE \
-dg $gene \
-o $gene

# Get into the parent directory 
#cd $mydir
