# This script adds chromosomes to .bed
# files that miss them and removes 
# chrY and chrM

# Get the file without 'chr' in chromosome name
input=$1

if [ $# -ne 1 ]; then
	echo "Please provide the right number of arguments"
	echo "USAGE: add_chr_to_bed.sh [.bed FILE WITHOUT CHR NAMES]"
	exit
fi

# Add the chromosome names and remove chrY and chrM
cat $input | sed 's/^/chr/' | grep -vE 'chrY|chrM' > ${input%.bed}.chr.bed
