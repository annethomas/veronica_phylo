#!/bin/bash

main_dir=~/HybPiper/genbank_compare
gene_dir="$main_dir"/alignments/final

out_dir="$main_dir"/alignments/trimmed
mkdir $out_dir

cd $gene_dir

for file in *.fasta 
do 
filebase=${file%.fasta}
echo $filebase
trimal -gt 0.3 -in $file -out ${out_dir}/${filebase}_trm.fasta
done 

