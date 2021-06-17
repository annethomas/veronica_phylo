#!/bin/bash
compare_dir=~/HybPiper/genbank_compare
tree_dir="$compare_dir"/iqtree

mkdir $tree_dir

cd $tree_dir


/scripts/csmit "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${tree_dir}/Veronica_genbank_concat5.phy -spp ${tree_dir}/Veronica_genbank_concat5_gene_pos.nex -m MF+MERGE" 
#Veronica_genbank_concat5_gene_pos.nex.best_scheme

#/scripts/csmit -b -c 16 -m 8G "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${main_dir}/${exon_filebasename}.phy -spp ${iq_dir}/${exon_filebasename}_gene_pos_iq.nex -m MFP+MERGE -rcluster 10