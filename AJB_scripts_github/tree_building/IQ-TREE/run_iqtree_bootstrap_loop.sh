#!/bin/bash

iq_dir=~/HybPiper/genbank_compare/iqtree

filebasename=$1       ##e.g. Veronica_gebank_concat5
base_dir="$iq_dir"/$1

num_bootstraps_start=$2 ##specify starting point for number of sets of 10 (eg 2 will do 2:total num of sets)
num_bootstraps_end=$3 ##specify end number of sets of 10

cd $base_dir

for i in $(seq $num_bootstraps_start $num_bootstraps_end)
do echo "submitting bootstrap set $i"
/scripts/csmit -b -c 8 "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${base_dir}/${filebasename}.phy -spp ${base_dir}/${filebasename}*.best_scheme.nex -b 10 -nt AUTO -ntmax 8 -pre ${filebasename}_bs$i"
echo "~/HybPiper/scripts/iqtree-1.6.12-Linux/bin/iqtree -s ${base_dir}/${filebasename}.phy -spp ${base_dir}/${filebasename}*.best_scheme.nex -b 10 -nt AUTO -ntmax 8 -pre ${filebasename}_bs$i"
sleep 1s
done

