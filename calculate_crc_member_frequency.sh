#!/bin/bash

# Use the output of CRCmapper as input to 
input=$1

no_networks=$(wc -l $input | cut -d' ' -f1)
echo $no_networks

cut -f1 $input | sed 's/\[//' | sed 's/\]//' | sed "s/'//g" | tr '\n' ',' | tr ',' '\n' | sed 's/^[[:space:]]*//' | sort | uniq -c | sort -k1,1nr | awk -v no_net=$no_networks '{print $2"\t"$1"\t"no_net}' > CRC_member_frequencies.tsv
