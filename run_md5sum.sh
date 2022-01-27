#!/bin/bash

file=$1

sum=$(md5sum $file)
echo $sum >> md5sums.txt
