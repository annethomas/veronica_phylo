#!/bin/bash

main_dir=~/HybPiper/genbank_compare
FNA_dir=${main_dir}/sequences
out_dir=${main_dir}/alignments/final

mkdir $out_dir

cd $FNA_dir

for file in *.fasta
   do 
   echo $file
   filebase=${file%.fasta}
   echo $filebase
   mafft --localpair --maxiterate 1000 $file > ${out_dir}/"$filebase"_aln.fasta
done 

mafft --adjustdirection ~/HybPiper/genbank_compare/sequences/Veronica_genbank_rpoB_trnC.fasta > ~/HybPiper/genbank_compare/alignments/final/Veronica_genbank_rpoB_trnC_aln.fasta
sed -i 's/_R_//g' ~/HybPiper/genbank_compare/alignments/final/Veronica_genbank_rpoB_trnC_aln.fasta
