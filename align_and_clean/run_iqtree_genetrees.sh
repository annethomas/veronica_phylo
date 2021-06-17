#!/bin/bash

main_dir=/production/at820/genbank_subset
compare_dir=~/HybPiper/genbank_compare

filebase=$1 ## e.g. exons_inframe_gsubset_trm or the cleaned version preferably

align_dir=${main_dir}/$1
tree_dir=${compare_dir}/iqtree/genetrees/$1

subset_file=$2 ## input gene subset file as argument (full path)
echo $subset_file

suffix=$3

mkdir -p $tree_dir

cd $align_dir

while read gene
do 
echo $gene
~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${align_dir}/${gene}"$suffix" -nt AUTO -ntmax 6 -pre "$tree_dir"/$gene
done < $subset_file

#~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s /production/at820/genbank_subset/exons_inframe_gsubset_trm_cleaned/4471_gsubset_trm.phy -nt AUTO -ntmax 6 -pre 4471_test