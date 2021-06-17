#!/bin/bash

compare_dir=~/HybPiper/genbank_compare/
main_dir="$compare_dir"/genbank_subset
iq_dir=~/HybPiper/genbank_compare/iqtree

concatname=$1 ## e.g. Veronica_229supercontigs_gsubset_concat
treealndir=$2 ## e.g. supercontigs_gsubset_trm_cleaned_tiptrim
basedir=$3 ## iqtree or astral
##optional
filter_file=$4 ## e.g. supercontigs_gsubset_trm_cleaned_tiptrim_secondpass.txt
type=$5 # e.g. secondpass
in_tree_suffix=$6
in_aln_suffix=$7

if [ $basedir == astral ] 
then concatfile="$compare_dir"/"$basedir"/"$concatname"/"$concatname".astral.tre
else concatfile="$compare_dir"/"$basedir"/"$concatname"/"$concatname"_allbstrees.suptree
fi

genetree_dir="$iq_dir"/genetrees/$treealndir
genealn_dir="$main_dir"/$treealndir
stat_dir=~/HybPiper/genbank_compare/stats/concordance_factors/$basedir
mkdir $stat_dir

## if the filter file is supplied, create a new folder with filtered genes
if [ ! -z "$filter_file" ]
  then
    filter_file=$main_dir/genelists/$filter_file
    echo "filtering trees"
    ~/HybPiper/scripts/utils/filter_gene_files.sh $genetree_dir "$genetree_dir"_"$type" $in_tree_suffix $filter_file
    genetree_dir="$genetree_dir"_"$type"
    echo "new gene tree directory: " $genetree_dir

    echo "filtering alignments"
    ~/HybPiper/scripts/utils/filter_gene_files.sh $genealn_dir "$genealn_dir"_"$type" $in_aln_suffix $filter_file
    genealn_dir="$genealn_dir"_"$type"
    echo "new aln directory: " $genealn_dir

    cd $genetree_dir

    #Rscript ~/HybPiper/genbank_compare/scripts/rename_tips.R $genetree_dir .treefile
    #Rscript ~/HybPiper/genbank_compare/scripts/rename_tips.R $genealn_dir .fasta

    treefilename="$genetree_dir"/"$concatname"_genetrees.treesfile
    cat *.treefile > $treefilename


else echo "Using supplied directories"
fi




treefilename="$genetree_dir"/"$concatname"_genetrees.treesfile
#cat *.treefile > $treefilename


cd $stat_dir
if [ $basedir == astral ] 
then ~/iqtree-2.0.6-Linux/bin/iqtree2 -t $concatfile --gcf $treefilename --prefix "$concatname"
else ~/iqtree-2.0.6-Linux/bin/iqtree2 -t $concatfile --gcf $treefilename -p $genealn_dir --scf 100 --prefix "$concatname"
fi
