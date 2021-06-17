#!/bin/bash

main_dir=/production/at820
compare_dir=~/HybPiper/genbank_compare

align_dir="$main_dir"/genbank_subset/$1 #(exon_gsubset_aln, etc)
stat_dir="$compare_dir"/stats/"$1"_stats
echo $stat_dir


## get the directory basename and remove any trailing slashes
out_prefix=$(echo "$1" | sed 's:/*$::') 
echo $out_prefix

mkdir $stat_dir
echo "gene,avg_identity" > $stat_dir/"$out_prefix"_avg_identity.csv
cd $align_dir


##get the generic alignment filename suffix without the gene
files=(*.fasta) #risks the suffix not including fasta
infile_suffix=${files[0]:4} 
echo $infile_suffix



while read gene
do 
#echo $gene
trimal -in ${align_dir}/"$gene""$infile_suffix" -sident > ${stat_dir}/${gene}.ident_stats.txt
avg_ident=$(grep -oP "(?<=^## AverageIdentity[[:blank:]])[0-9]\.[0-9]+" ${stat_dir}/${gene}.ident_stats.txt)
echo "$gene,$avg_ident"
echo "$gene,$avg_ident" >> $stat_dir/"$out_prefix"_avg_identity.csv
done < ${main_dir}/genelists/genes_firstpass_270.txt

