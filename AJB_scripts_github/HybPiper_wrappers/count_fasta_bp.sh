#!/bin/bash

# count total base pairs in all samples' output as collated by retrive_sequences.py (unaligned fasta for each gene or flanking introns))
cd /production/at820/supercontig
bp_total=0

# for file in *.FNA # for the exons
for file in *.fasta # for the introns
do
   bp=$(grep -v ">" $file | wc | awk '{print $3-$1}')
   let "bp_total = bp_total + bp"
done

echo $bp_total