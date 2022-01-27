#!/bin/bash

input=$1
feature=$2

cat $input | awk -v var=$feature '{print $1"\t"$2"\t"$3"\t"var"_"NR"\t"".""\t"$5}' | sortBed > ${input%.txt}.bed
