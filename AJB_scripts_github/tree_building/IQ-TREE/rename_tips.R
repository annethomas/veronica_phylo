## rename labels from sample ID to species name

hyb_dir="~/HybPiper/genbank_compare/"

in_dir=commandArgs(trailingOnly=TRUE)[1] # full filepath e.g. ...exons_inframe_gsubset_trm_cleaned/trimmed_trees 
suffix=commandArgs(trailingOnly=TRUE)[2] #e.g. .treefile or .fasta

#in_dir=file.path(hyb_dir,"iqtree/genetrees/",basename)
print(in_dir)

libs=read.csv(file.path(hyb_dir,"all_libs.csv"),stringsAsFactors = FALSE)
print(paste0(suffix,"$"))	
files=list.files(in_dir,pattern=paste0(suffix,"$"))
print(files)

setwd(in_dir)
if(grepl("fa",suffix)){
	print("reading alignments")
	
	for(file in files){
		print(file)
		aln=ape::read.dna(file,format="fasta")
		labs=data.frame(lab=dimnames(aln)[[1]])
		if(any(grepl("V",dimnames(aln)[[1]]))){
			print("labels already renamed to species")
			next
		}

		## remove gene numbers from some labels (introns and supercontigs; won't do anything otherwise)
		labs$lab=gsub("-[0-9]{4}","",labs$lab)

		## join to species-sampleid dataframe to match them up
		labs=dplyr::left_join(labs,libs,by=c("lab"="Sample"))

		## add underscores
		labs$Species=gsub(" ","_",labs$Species)

		## take out the extra id on the two doubled species (since we're only using one of each and genbank only has the species name)
		labs$Species=gsub("-.*[0-9]+$","",labs$Species)

		## replace old labels and write to file
		dimnames(aln)[[1]]=labs$Species
		#print(dimnames(aln)[[1]])
		ape::write.dna(aln,file,format="fasta",nbcol=-1,colsep="")
	}

} else if(grepl("tr",suffix)){
	print("reading trees")
		
	for(tree_file in files){
		print(tree_file)
		tree=ape::read.tree(tree_file)
		
		if(any(grepl("V",tree$tip.label))){
			stop("labels already renamed to species")
			next
		}

		tips=data.frame(tip=tree$tip.label)
		tips$tip=gsub("-[0-9]{4}","",tips$tip)
		tips=dplyr::left_join(tips,libs,by=c("tip"="Sample"))
		tips$Species=gsub(" ","_",tips$Species)
		tips$Species=gsub("-.*[0-9]+$","",tips$Species)
		tree$tip.label=tips$Species
		#print(tree$tip.label)
		ape::write.tree(tree,tree_file)
	}
} else{
	print("Can't tell what filetype this is")
}

#[ultimately build this in the the genetree generating script?]