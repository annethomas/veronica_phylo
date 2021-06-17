#!/bin/bash
compare_dir=~/HybPiper/genbank_compare
tree_dir="$compare_dir"/iqtree

aln_base=$1
aln_dir="$tree_dir"/$1 ## e.g. Veronica_198exons_gsubset_concat

cd $aln_dir

echo "working directory:"
pwd

#~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${aln_base}.phy -spp ${aln_base}_gene_pos.nex -m MF+MERGE -rcluster 10 -nt AUTO -ntmax 4 

/scripts/csmit -c 16 -m 8G "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${aln_base}.phy -spp ${aln_base}_gene_pos.nex -m MF+MERGE -rcluster 10 -nt AUTO -ntmax 16" 


#/scripts/csmit -b -c 16 -m 8G "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${main_dir}/${exon_filebasename}.phy -spp ${iq_dir}/${exon_filebasename}_gene_pos_iq.nex -m MFP+MERGE -rcluster 10