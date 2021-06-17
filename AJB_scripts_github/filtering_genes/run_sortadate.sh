#!/bin/bash

compare_dir=~/HybPiper/genbank_compare/
statdir="$compare_dir"/stats/sortadate

mkdir $statdir

ingenetreedir="$compare_dir"/iqtree/genetrees/$1 ##e.g. exons_gsubset_trm_cleaned_tiptrim
genetreedir="$ingenetreedir"_rooted
outfile_base=$1

ref_tree=$2 #e.g. Veronica_173exons_gsubset_concat

## root the trees
~/HybPiper/genbank_compare/scripts/post_tree_analysis/reroot_trees_phyx.sh $outfile_base $ref_tree
input_tree="$compare_dir"/iqtree/"$ref_tree"/"$ref_tree"_allbstrees.suptree.rooted 

cd ~/SortaDate-master

# Get the root-to-tip variance and tree length
python src/get_var_length.py "$genetreedir" --flend .treefile --outf "$statdir"/"$outfile_base"_variance.txt --outg Veronica_chamaedrys

# Get the bipartition support 

python src/get_bp_genetrees.py $genetreedir "$input_tree" --flend .treefile --outf "$statdir"/"$outfile_base"_bpsupport.txt

# Combine the results from these two runs 
python src/combine_results.py "$statdir"/"$outfile_base"_variance.txt "$statdir"/"$outfile_base"_bpsupport.txt --outf "$statdir"/"$outfile_base"_combined.txt

## save separate file with the NA rows removed
sed '/NA/d' "$statdir"/"$outfile_base"_combined.txt > "$statdir"/"$outfile_base"_combined_noNA.txt

# Sort and get the list of the good genes (assuming the order is root-to-tip_var, treelength, bipartition)
python src/get_good_genes.py "$statdir"/"$outfile_base"_combined_noNA.txt --max 270 --order 3,1,2 --outf "$statdir"/"$outfile_base"_sortedgenes.txt


## debugging
# pxrmt -t ~/HybPiper/genbank_compare/iqtree/genetrees/exons_gsubset_trm_cleaned_tiptrim/4471.treefile -n Veronica_chamaedrys | pxlstr
# tree #: 0
# rooted: false
# binary: true
# nterminal: 77
# ninternal: 75
# branch lengths: true
# rttipvar: NA
# treelength: 0.608264
# ultrametric: false
# rootheight: NA

## doesnt' contain chamaedrys (but does contain perfoliata)
#pxrr -t ~/HybPiper/genbank_compare/iqtree/genetrees/exons_gsubset_trm_cleaned_tiptrim/4942.treefile -r -g Veronica_chamaedrys,Veronica_perfoliata -o testroot_4942.tre
## contains both
#pxrr -t ~/HybPiper/genbank_compare/iqtree/genetrees/exons_gsubset_trm_cleaned_tiptrim/4471.treefile -r -g Veronica_chamaedrys,Veronica_perfoliata -o testroot3_4471.tre

#grep -L 'Veronica_chamaedrys' *.treefile
#grep -L 'Veronica_chamaedrys\|Veronica_perfoliata' *.treefile