#!/bin/bash

main_dir=~/HybPiper/
in_dir=$1 # full file path e.g. ~/HybPiper/alignments/supercontigs_gt70_cleaned_tiptrim
out_dir=$2 ## e.g. ~/HybPiper/alignments/supercontigs_tiptrim_sortadateRank20
in_suffix=$3 ## e.g. _supercontig_trimmed_gt70.tt.fasta
filter_file=$4 ## e.g. ~/HybPiper/genelists/sortadate_wtrank_20genes.txt

echo $filter_file
dos2unix $filter_file

mkdir $out_dir

while read gene
do 
echo $gene
cp ${in_dir}/"$gene""$in_suffix" $out_dir
done < $filter_file

## ~/HybPiper/scripts/utils/filter_gene_files.sh ~/HybPiper/alignments/supercontigs_gt70_cleaned_tiptrim ~/HybPiper/alignments/supercontigs_tiptrim_sortadateRank20 _supercontig_trimmed_gt70.tt.fasta ~/HybPiper/genelists/sortadate_wtrank_20genes.txt
## ~/HybPiper/scripts/utils/filter_gene_files.sh ~/HybPiper/alignments/supercontigs_gt70_cleaned_tiptrim ~/HybPiper/alignments/supercontigs_tiptrim_sortadateRank20 _supercontig_trimmed_gt70.tt.fasta ~/HybPiper/genelists/sortadate_wtrank_20genes.txt
