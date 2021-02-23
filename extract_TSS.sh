# Extract the TSS from the refGene file
input=$1

# Convert the refGene file to .bed
#awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,0,$5}' $input > ucsc_canonical.bed

# Extract the TSS (depends on the strand)
awk 'BEGIN{FS=OFS="\t"}($6=="+"){print $1,$2,$2,$4,$5,$6}' $input > ${input%.bed}_TSS.bed
awk 'BEGIN{FS=OFS="\t"}($6=="-"){print $1,$3,$3,$4,$5,$6}' $input >> ${input%.bed}_TSS.bed
