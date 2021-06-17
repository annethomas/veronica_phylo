#!/bin/bash

main_dir=/production/at820

cd "$main_dir"

python /applications/HybPiper/get_seq_lengths.py "$main_dir"/Angiosperms353_targetSequences.fasta "$main_dir"/veronica_seq_prefixes.txt dna > "$main_dir"/stats/gene_lengths.txt

python /applications/HybPiper/hybpiper_stats.py "$main_dir"/stats/gene_lengths_20200311.txt "$main_dir"/veronica_seq_prefixes.txt > "$main_dir"/stats/hybpiper_stats.txt
