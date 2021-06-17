#!/bin/bash
## divide the overall genelist into subsets of specified size and pass in file to run_iqtree_genetrees to run as a job
#~/HybPiper/genbank_compare/scripts/run_iqtree_genetrees_loop_wrapper.sh exons_gsubset_trm_cleaned_tiptrim /production/at820/genbank_subset/exons_gsubset_trm_cleaned_tiptrim_secondpass3.txt _gsubset_trm.tt.fasta 20
##~/HybPiper/genbank_compare/scripts/run_iqtree_genetrees_loop_wrapper.sh introns_gsubset_trm_cleaned_tiptrim /production/at820/genbank_subset/introns_gsubset_trm_cleaned_tiptrim_secondpass.txt _gsubset_trm.tt.fasta 20

main_dir=/production/at820/genbank_subset
compare_dir=~/HybPiper/genbank_compare

filebase=$1 ## e.g. exons_inframe_gsubset_trm or the cleaned version preferably
subset_file=$2 ## overall input gene subset file as argument (full path)
suffix=$3 ## e.g. _gsubset_trm.tt.fasta
subset_num=$4 ## number of files to include in one subset



nlines=$(wc -l < $subset_file)
let "roundup = $subset_num - 1"
let "niter = ( $nlines + $roundup ) / $subset_num"

for i in $(seq 1 $niter)
	do echo $i
	let "start = ( $i * $subset_num ) - ( $subset_num - 1 )"
	let "end = $i * $subset_num"
	subset_fn="$subset_file"_subset"$i".txt
	
	echo $subset_fn
	echo "$start to $end"
	
	sed -n "$start","$end"p $subset_file > $subset_fn
	
	/scripts/csmit -b -c 6 "~/HybPiper/genbank_compare/scripts/run_iqtree_genetrees.sh '$filebase' '$subset_fn' '$suffix'"
	sleep 1s
done




















