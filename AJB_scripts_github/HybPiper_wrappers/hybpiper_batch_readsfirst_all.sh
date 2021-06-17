#!/bin/bash
#fastqdir=~/HybPiper/FASTQ_Generation_2019-05-16/
fastqdir=~/HybPiper/FASTQ_2019-10/

# export PATH=$PATH:/applications/spades/bin:/applications/spades/share/spades:/applications/fastqc/fastqc_v0.11.5:/applications/bwa/bwa-0.7.17/bwa:/usr/bin/exonerate:/usr/bin/parallel:/usr/bin/samtools/:/usr/bin/makeblastdb:/usr/bin/blastx:

echo $PATH
python /applications/HybPiper/reads_first.py --check

cd ${fastqdir}/fastq_trimmed_pairs/


while read name;
do /scripts/csmit -b -c 8 "/applications/HybPiper/reads_first.py -b ~/HybPiper/scripts/Angiosperms353_targetSequences.fasta -r ${name}_trimmed_R*.fastq --prefix /production/at820/${name} --bwa"
sleep 1s
done < ${fastqdir}/seq_prefixes.txt


