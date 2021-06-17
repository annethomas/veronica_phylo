#! /bin/bash

compare_dir=~/HybPiper/genbank_compare
data_dir=/production/at820/genbank_subset

input_treedir="$compare_dir"/iqtree/genetrees/$1
output_treedir="$compare_dir"/iqtree/genetrees/$2

input_alndir="$data_dir"/$1
output_alndir="$data_dir"/$2

filter_file=$data_dir/genelists/$3
dos2unix $filter_file

mkdir $output_treedir
mkdir $output_alndir

while read gene
do
	echo $gene
	cp "$input_treedir"/"$gene"*.treefile "$output_treedir"
	cp "$input_alndir"/"$gene"*.fasta "$output_alndir"
done < $filter_file