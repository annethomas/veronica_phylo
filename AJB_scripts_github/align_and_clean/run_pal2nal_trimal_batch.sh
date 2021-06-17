#!/bin/bash

## edit directory names if necessary
## add /applications/trimal/source to your PATH variable

main_dir=/production/at820/genbank_subset
FAA_dir=${main_dir}/exons_FAA_gsubset_aln/
FNA_dir=${main_dir}/exons_FNA_gsubset_aln/
inframe_dir=${main_dir}/exons_inframe_gsubset_aln/
trim_dir=${main_dir}/exons_inframe_gsubset_trm/

cd $main_dir
mkdir $inframe_dir
mkdir $trim_dir



while read gene
do
echo $gene
~/HybPiper/scripts/pal2nal.v14/pal2nal.pl -output fasta ${FAA_dir}/${gene}_gsubset_aln.fasta ${FNA_dir}/${gene}_gsubset_aln.fasta > ${inframe_dir}/${gene}_gsubset_infr.fasta
# these trimal settings will remove columns in the alignment with >50% gaps or similarity score of < .001, but won't remove more than 40%; it also removes sequences with < 50% good sites (good sites determined by resoverlap)
trimal -gt 0.5 -st .001 -cons 0.4 -resoverlap 0.5 -seqoverlap 50 -in ${inframe_dir}/${gene}_gsubset_infr.fasta -out ${trim_dir}/${gene}_gsubset_trm.fasta
done < /production/at820/genelists/genes_firstpass_270.txt
