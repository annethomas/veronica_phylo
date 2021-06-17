#!/bin/bash

main_dir=~/HybPiper/genbank_compare
seq_dir=${main_dir}/sequences

for file in *.fasta
do
	out_file=${file%.*}.csv

	grep "^>" $file | sed -e 's/>//' -e 's/_/,/' > $out_file
done
