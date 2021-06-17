#!/bin/bash

main_dir=/production/at820

cd $main_dir

while read name;
do
echo $name
python /applications/HybPiper/paralog_investigator.py $name 2>> "$main_dir"/stats/paralog_warnings.txt
done < veronica_seq_prefixes.txt
