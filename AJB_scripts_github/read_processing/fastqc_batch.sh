fastq_dir=~/HybPiper/FASTQ_2019-10/
fastq_dir_version=${fastq_dir}/fastq_concat_pretrim ## change to run for trimmed fastq files
fastqc_dir=${fastq_dir}/fastqc_pretrim 


## run fastqc on all reads
cd $fastq_dir_version
fastqc *.fastq.gz

## move fastqc files to a different dir (I’ve had problems with the –outdir option in fastqc)
mkdir $fastqc_dir
mv ${fastq_dir_version}/*fastqc* $fastqc_dir

## (output includes html summaries)
