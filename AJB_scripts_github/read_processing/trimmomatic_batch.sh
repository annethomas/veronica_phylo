#!/bin/bash

#fastqdir=~/HybPiper/FASTQ_Generation_2019-05-16/
fastqdir=~/HybPiper/FASTQ_2019-10/
mkdir ${fastqdir}/fastq_trimmed

cd ${fastqdir}/fastq_concat

gunzip *.gz

while read prefix;
do 
echo $prefix
java -jar /applications/trimmomatic/Trimmomatic-0.34/trimmomatic-0.34.jar PE ${prefix}_R1.fastq ${prefix}_R2.fastq -baseout ${fastqdir}/fastq_trimmed/${prefix}_trimmed.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
done < ${fastqdir}/seq_prefixes.txt

## move paired files to separate directory
mkdir ${fastqdir}/fastq_trimmed_pairs
mv ${fastqdir}/fastq_trimmed/*P.fastq ${fastqdir}/fastq_trimmed_pairs


