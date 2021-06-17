#!/bin/bash

#fastqdir=~/HybPiper/FASTQ_Generation_2019-05-16/
fastqdir=~/HybPiper/FASTQ_2019-10/

mkdir ${fastqdir}/fastq_concat
cd ${fastqdir}/fastq_pre_concat

pwd

while read prefix;
do
echo $prefix
reads1=$(grep ${prefix}_.*_R1 ${fastqdir}/sequence_filenames.txt)
cat $reads1 > ${fastqdir}/fastq_concat/${prefix}_R1.fastq.gz
reads2=$(grep ${prefix}_.*_R2 ${fastqdir}/sequence_filenames.txt)
cat $reads2 > ${fastqdir}/fastq_concat/${prefix}_R2.fastq.gz
done < ${fastqdir}/seq_prefixes.txt
