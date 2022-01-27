#!/bin/bash

input=$1

cat $input | sort | uniq | sortBed > ${input%.bed}.uniq.sorted.bed
