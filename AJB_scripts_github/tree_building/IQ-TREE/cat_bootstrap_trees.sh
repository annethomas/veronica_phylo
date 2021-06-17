#!/bin/bash

## first run run_MLforBS.sh as job

main_dir=/production/at820/genbank_subset
iq_dir=~/HybPiper/genbank_compare/iqtree

filebasename=$1 ## e.g. Veronica_227supercontigs_gsubset_concat
tree_dir="$iq_dir"/"$filebasename"


cd ${tree_dir}
cat *.boottrees > ${filebasename}_allbstrees
~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -con -t ${filebasename}_allbstrees

## replaced by run_MLforBS.sh so that it can be run before bootstraps are done
#~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -nt AUTO -ntmax 12 -s ${tree_dir}/${filebasename}.phy -spp ${tree_dir}/${filebasename}*.best_scheme.nex -pre ${filebasename}_MLforBS

~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -sup ${filebasename}_MLforBS.treefile -t ${filebasename}_allbstrees


