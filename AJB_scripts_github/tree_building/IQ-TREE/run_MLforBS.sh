#!/bin/bash

iq_dir=~/HybPiper/genbank_compare/iqtree

filebasename=$1 ## e.g. Veronica_227supercontigs_gsubset_concat
tree_dir="$iq_dir"/"$filebasename"


cd ${tree_dir}

/scripts/csmit -c 12 "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -nt AUTO -ntmax 12 -s ${tree_dir}/${filebasename}.phy -spp ${tree_dir}/${filebasename}*.best_scheme.nex -pre ${filebasename}_MLforBS"


