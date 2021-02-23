# Get the input
input=$1

cat $input | awk -F"\t" '{print $2"\t"$4"\t"$5"\t"$1"\t""0""\t"$3}' | sortBed > ${input%.txt}.bed
