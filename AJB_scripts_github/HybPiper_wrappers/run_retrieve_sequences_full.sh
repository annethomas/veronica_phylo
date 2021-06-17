#!/bin/bash

main_dir=/production/at820/

mkdir ${main_dir}/exons_FNA/
cd ${main_dir}/exons_FNA/
pwd
/scripts/csmit -b "python /applications/HybPiper/retrieve_sequences.py ~/HybPiper/scripts/Angiosperms353_targetSequences.fasta $main_dir dna"
sleep 1s

mkdir ${main_dir}/exons_FAA/
cd ${main_dir}/exons_FAA/
pwd
/scripts/csmit -b "python /applications/HybPiper/retrieve_sequences.py ~/HybPiper/scripts/Angiosperms353_targetSequences.fasta $main_dir aa"
sleep 1s

## for amino acids
while read line
do
sed -i 's/*/X/g' $line
done < <(ls)


mkdir ${main_dir}/introns/
cd ${main_dir}/introns/
pwd
/scripts/csmit -b "python /applications/HybPiper/retrieve_sequences.py ~/HybPiper/scripts/Angiosperms353_targetSequences.fasta $main_dir intron"
sleep 1s

mkdir ${main_dir}/supercontigs/
cd ${main_dir}/supercontigs/
pwd
/scripts/csmit -b "python /applications/HybPiper/retrieve_sequences.py ~/HybPiper/scripts/Angiosperms353_targetSequences.fasta $main_dir supercontig"
