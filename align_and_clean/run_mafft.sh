#!/bin/bash

main_dir=/production/at820/
gene_dir="$main_dir"/genbank_subset/$1_genbanksubset #(exons_FNA, exons_FAA, introns, supercontig)

out_dir="$gene_dir"_aln/
mkdir $out_dir

cd $gene_dir

while read gene
do 
echo $gene
mafft --auto ${gene}_gsubset.fasta > ${out_dir}/${gene}_gsubset_aln.fasta
done < ${main_dir}/genelists/genes_firstpass_270.txt
