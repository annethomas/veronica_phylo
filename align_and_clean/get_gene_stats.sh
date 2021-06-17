#!/bin/bash

main_dir=/production/at820
compare_dir=~/HybPiper/genbank_compare

subset_dir="$main_dir"/genbank_subset/$1 #(exons_FNA_gsubset, exons_FAA_gsubset, introns_gsubset, supercontig_gsubset, exon_gsubset_aln, etc)



cd $subset_dir

for file in *.fasta
	do 
	if [ $(grep " " $file | wc -l) -gt 0 ]
	then 
		Rscript ~/HybPiper/genbank_compare/scripts/rewrite_aln_nocol.R $file
		echo "rewrote "$file" to remove spaces"
	fi

done

seqkit stat --tabular --all *.fasta > "$compare_dir"/stats/"$1"_stats/"$1"_stats.txt



#main_dir=/production/at820/
#seq_dir=${main_dir}/$1

#cd $seq_dir

#while read gene
#do 
#echo $gene
#mafft --localpair --maxiterate 1000 ${gene}.FNA > ${out_dir}${gene}_aln.FNA
#done < ${main_dir}/genes_firstpass_270.txt

#if [ $(grep " " $file | wc -l) -gt 0 ]; then Rscript ~/HybPiper/genbank_compare/scripts/rewrite_aln_nocol.R $file; echo "rewrote "$file" to remove spaces"; fi
	
