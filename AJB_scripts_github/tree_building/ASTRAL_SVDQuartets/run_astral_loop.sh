#!/bin/bash

mode=$1 # leave blank to do both; gcf if only gcf

if [ -z $mode ]
then
echo 1
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_229supercontigs_gsubset_concat supercontigs_gsubset_trm_cleaned_tiptrim
echo 2
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_167introns_gsubset_concat introns_gsubset_trm_cleaned_tiptrim
echo 3
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_173exons_gsubset_concat exons_gsubset_trm_cleaned_tiptrim #exons_gsubset_trm_cleaned_tiptrim_secondpass.txt .treefile

echo 4
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_supercontigs_intersect_concat supercontigs_intersect
echo 5
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_introns_intersect_concat introns_intersect
echo 6
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_exons_intersect_concat exons_intersect

 echo 7
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_fix150supercontigs_gsubset_concat supercontigs_fix150
 echo 8
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_fix116introns_gsubset_concat introns_fix116
 echo 9
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_fix96exons_gsubset_concat exons_fix96

echo 10
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_med80supercontigs_gsubset_concat supercontigs_med80
echo 11
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_med63introns_gsubset_concat introns_med63
echo 12
~/HybPiper/genbank_compare/scripts/run_astral.sh Veronica_med51exons_gsubset_concat exons_med51
fi

## gcf
echo "running gcf"

#echo 1
#~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_229supercontigs_gsubset_concat supercontigs_gsubset_trm_cleaned_tiptrim astral
#echo 2
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_167introns_gsubset_concat introns_gsubset_trm_cleaned_tiptrim astral introns_gsubset_trm_cleaned_tiptrim_secondpass.txt secondpass .treefile _gsubset_trm.tt.fasta
echo 3
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_173exons_gsubset_concat exons_gsubset_trm_cleaned_tiptrim astral exons_gsubset_trm_cleaned_tiptrim_secondpass.txt secondpass .treefile _gsubset_trm.tt.fasta

echo 4
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_supercontigs_intersect_concat supercontigs_intersect astral
echo 5
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_introns_intersect_concat introns_intersect astral
echo 6
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_exons_intersect_concat exons_intersect astral

 echo 7
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_fix150supercontigs_gsubset_concat supercontigs_fix150 astral
 echo 8
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_fix116introns_gsubset_concat introns_fix116 astral
 echo 9
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_fix96exons_gsubset_concat exons_fix96 astral

echo 10
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_med80supercontigs_gsubset_concat supercontigs_med80 astral
echo 11
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_med63introns_gsubset_concat introns_med63 astral
echo 12
~/HybPiper/genbank_compare/scripts/post_tree_analysis/run_iqtree_gcf.sh Veronica_med51exons_gsubset_concat exons_med51 astral


