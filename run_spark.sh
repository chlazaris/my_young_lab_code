#!/bin/bash

# Get the current directory
mydir=$pwd

# Get into the other dir
cd input

# Get the required inputs
coords=$1
gene=$2

#Code used to generate this plot:
python /lab/solexa_young/lazaris/tools/SparK/SparK.py \
-pr $coords \
-cf H3K27ac.bdg MED1.bdg BRD4.bdg CO_IKZF1.bdg IKZF1.bdg \
-gtf gencode.v19.annotation.gtf \
-gl H3K27ac MED1 BRD4 CO_IKZF1 IKZF1 \
-gs yes \
-sm 20 \
-cs 10 10 10 10 10 \
-bed MM1S_superenhancer_peaks.bed IKZF1_peaks.bed \
-f 0000FF \
-bedcol 000000 000000 \
-bedlab super-enhancer IKZF1 \
-dg $gene \
-o ../results/$gene

# Get into the parent directory 
cd $mydir
