#!/bin/bash

# Get the current date
now=$(date +"%m-%d-%y")

# Generate tableau input
grep . */used_disk_space_in_GB.txt | sed 's/\//	/' | sed 's/\:/	/' | cut -d'	' -f1,3 > ${now}_disk_space_for_tableau.tsv 
