#!/bin/bash

main_dir=~/HybPiper/genbank_compare/
concatname=$1 	## e.g. Veronica_229supercontigs_gsubset
genetree_dir=${main_dir}/iqtree/genetrees/$2 	## e.g. supercontigs_gsubset_trm_cleaned_tiptrim
## optional:
filter_file=$3 ## relative filepath, e.g. supercontigs_gsubset_trm_cleaned_tiptrim_secondpass.txt
filter_suffix=$4 ## e.g. .treefile

astral_dir=${main_dir}/astral/$concatname
mkdir $astral_dir

treesfilename="$genetree_dir"/"$concatname"_genetrees.treesfile
cd $genetree_dir

## if the filter file is supplied, make concatenated genetree file; otherwise assume it already exists in genetree_dir
if [ ! -z "$filter_file" ]
then
	filter_file="$main_dir"/genbank_subset/genelists/"$filter_file"
	echo "$filter_file"
	if [ -f "$treesfilename" ]
	then
		echo "concatenated genetree file already exists"
		exit 1
	else
		echo "creating $treesfilename from $filter_file"
		dos2unix "$filter_file"
		while read gene
		do 
			echo $gene
			cat "$gene""$filter_suffix" >> $treesfilename
		done < $filter_file
	fi
fi


cd $astral_dir

java -jar ~/Astral/astral.5.6.3.jar -i $treesfilename -o $concatname.astral.tre 2> $concatname.astral.log



