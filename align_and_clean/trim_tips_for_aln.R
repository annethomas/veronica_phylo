#Rscript ~/HybPiper/genbank_compare/scripts/trim_tips_for_aln.R ~/HybPiper/genbank_compare/iqtree/genetrees/exons_inframe_gsubset_trm_cleaned/trimmed_trees .treefile /production/at820/genbank_subset/exons_inframe_gsubset_trm_cleaned .fasta /production/at820/genbank_subset/exons_gsubset_trm_cleaned_tiptrim S4D H2L
#Rscript ~/HybPiper/genbank_compare/scripts/trim_tips_for_aln.R ~/HybPiper/genbank_compare/iqtree/genetrees/supercontigs_gsubset_trm_cleaned/trimmed_trees .treefile /production/at820/genbank_subset/supercontigs_gsubset_trm_cleaned .fasta /production/at820/genbank_subset/supercontigs_gsubset_trm_cleaned_tiptrim S4D H2L
#Rscript ~/HybPiper/genbank_compare/scripts/trim_tips_for_aln.R ~/HybPiper/genbank_compare/iqtree/genetrees/supercontigs_gsubset_trm_cleaned/trimmed_trees/test .treefile /production/at820/genbank_subset/supercontigs_gsubset_trm_cleaned .fasta ~/HybPiper/genbank_compare/iqtree/genetrees/supercontigs_gsubset_trm_cleaned/trimmed_trees/test S4D H2L
#Rscript ~/HybPiper/genbank_compare/scripts/trim_tips_for_aln.R ~/HybPiper/genbank_compare/iqtree/genetrees/introns_gsubset_trm_cleaned/trimmed_trees/ .treefile /production/at820/genbank_subset/introns_gsubset_trm_cleaned .fasta /production/at820/genbank_subset/introns_gsubset_trm_cleaned_tiptrim S4D H2L


library(ape)

args=commandArgs(trailingOnly=TRUE)
trim_tree_dir=args[1]
tree_suffix=args[2]
aln_dir=args[3]
aln_suffix=args[4]
out_aln_dir=args[5]
outgroups=args[6:length(args)]

dir.create(out_aln_dir)

tree_files=list.files(trim_tree_dir,pattern=paste0(tree_suffix,"$"))
print(tree_files)
print(outgroups)
aln_files=list.files(aln_dir,pattern=aln_suffix)
stats_df=data.frame(gene=character(0),num_seqs_removed=numeric(0),seqs=character(0),outgroups_rescued=character(0),stringsAsFactors=FALSE)

for(file in tree_files){
	print(file)
	
	tree1_fn=file.path(trim_tree_dir,file)
	tree2_fn=file.path(trim_tree_dir,paste0(file,".tt"))

	if(!file.exists(tree1_fn)){
		print(paste(tree1_fn,"is missing"))
		next
	} else if(!file.exists(tree2_fn)){
		print(paste(tree2_fn,"is missing"))
		stats_row=data.frame(gene=gene,num_seqs_removed="?",seqs="",outgroups_rescued="",stringsAsFactors=FALSE)
	
		print(stats_row)
		stats_df=rbind(stats_df, stats_row)

		next
	} 


	gene=substr(file,1,4)
	
	tree1=ape::read.tree(file.path(trim_tree_dir,file))
	tree1.tips=tree1$tip.label
	#print(tree1.tips)

	tree2=ape::read.tree(file.path(trim_tree_dir,paste0(file,".tt")))
	tree2.tips=tree2$tip.label
	#print(tree2.tips)

	seqs_to_remove=setdiff(tree1.tips,tree2.tips)
	outgroups_rescued = c()
	if(length(outgroups)>0){
		print("checking outgroups")
		if(grepl("supercontig", aln_dir) | grepl("intron",aln_dir)){
			print("fixing outgroups")
			outgroups_gene=paste(outgroups,gene,sep="-")
			print(outgroups_gene)
		}
		outgroups_rescued = outgroups_gene[which(outgroups_gene %in% seqs_to_remove)]
		print(outgroups_rescued)
		if(length(outgroups_rescued)>0){
			seqs_to_remove=seqs_to_remove[-which(seqs_to_remove %in% outgroups_gene)]
			print(paste("rescued",outgroups_rescued))
		}
	}
	
	print(seqs_to_remove)

	aln_file=grep(gene,aln_files,value=TRUE)

	if(length(aln_file)>1){
		stop(paste("more than one alignment for",gene)) 
	}
	newaln_fn=file.path(out_aln_dir,paste0(substr(aln_file,1,nchar(aln_file)-nchar(aln_suffix)),".tt",aln_suffix))

	if(length(seqs_to_remove)==0){
		print("no seqs to remove, copying original alignment")
		file.copy(file.path(aln_dir,aln_file),newaln_fn)
	} else{
		aln=ape::read.dna(file.path(aln_dir,aln_file),format="fasta")

		newaln=aln[-which(dimnames(aln)[[1]] %in% seqs_to_remove),]
		
		print(paste("new alignment file: ",newaln_fn))
		ape::write.dna(newaln,file=newaln_fn,format="fasta")
	}

	
	num_seqs_removed=length(seqs_to_remove)
	stats_row=data.frame(gene=gene,num_seqs_removed=num_seqs_removed,seqs=paste(seqs_to_remove,collapse=","),outgroups_rescued=paste(outgroups_rescued,collapse=","),stringsAsFactors=FALSE)
	
	print(stats_row)
	stats_df=rbind(stats_df, stats_row)

	#print("done")
}

stats_fn=file.path(out_aln_dir,"tiptrim_stats.txt")
write.table(stats_df,file=stats_fn,sep="\t",row.names=FALSE)
print(paste("stats written to", stats_fn))


