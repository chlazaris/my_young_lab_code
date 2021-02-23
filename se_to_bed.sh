#!/bin/bash

input=$1

cat $input | awk '$13==1' | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"0"\t""."}' | sortBed > ${input%.txt}.bed
