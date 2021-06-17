#!/bin/bash
fastqdir=~/HybPiper/FASTQ_Generation_2019-05-16/


echo $PATH
# python /applications/HybPiper/reads_first.py --check

cd /production/at820/


while read name;
do python /applications/HybPiper/intronerate.py --prefix ${name}
done < veronica_seq_prefixes.txt
