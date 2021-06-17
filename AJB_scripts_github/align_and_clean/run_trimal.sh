#!/bin/bash

main_dir=~/HybPiper/genbank_compare/
gene_dir="$main_dir"/genbank_subset/$1_gsubset_aln #(introns, supercontig)

out_dir="$main_dir"/genbank_subset/$1_gsubset_trm/
mkdir $out_dir

cd $gene_dir

while read gene
do
echo $gene
trimal -gt 0.7 -in ${gene}_gsubset_aln.fasta -out ${out_dir}/${gene}_gsubset_trm.fasta
done < ${main_dir}/genelists/genes_firstpass_270.txt
