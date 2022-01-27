input=$1

cat $input | \
	grep -v '^#' | awk '$3 == "gene"' | \
	grep gene_type=protein_coding | \
	grep -Ev 'chrM|chrY' | tr ';' ' ' | \
	awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$6"\t"$7}' | \
	sed 's/gene_id=//' | sed 's/gene_name=//' |  sed 's/\./	/' | \
	sed 's/\./ /' | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$6"\t"".""\t"$NF}' | \
	sort | uniq | sortBed  
