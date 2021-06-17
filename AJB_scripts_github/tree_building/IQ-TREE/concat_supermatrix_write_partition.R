## example:
## Rscript ~/HybPiper/genbank_compare/scripts/concat_supermatrix_write_partition.R exons_gsubset_trm_cleaned_tiptrim Veronica_173exons_gsubset_concat exons_gsubset_trm_cleaned_tiptrim_secondpass3.txt

if(!require(ape)) install.packages("ape",repos="https://cran.ma.imperial.ac.uk/")
if(!require(ape)) install.packages("dplyr",repos="https://cran.ma.imperial.ac.uk/")
library(ape)
library(dplyr)

## arguments: alignment directory name, output concat file basename, filter filename (eg supercontigs_secondpass.txt)

compare_dir= "~/HybPiper/genbank_compare/"
data_dir = file.path(compare_dir,"genbank_subset")

input_dir=commandArgs(trailingOnly=TRUE)[1] ## e.g. supercontigs_gsubset_trm_cleaned_tiptrim
outfile_basename=commandArgs(trailingOnly=TRUE)[2] ## e.g. Veronica_gsubset_195supercontigs_noparalogs
filter_file=commandArgs(trailingOnly=TRUE)[3] ## e.g. supercontigs_gsubset_trm_cleaned_tiptrim_secondpass


out_dir=file.path(compare_dir,"iqtree",outfile_basename)
dir.create(out_dir)

align_dir = file.path(data_dir,input_dir)
print(align_dir)
print(outfile_basename)
print(filter_file)

aln_suffix="*trm*"
aln_format="fasta"

libs=read.csv(file.path(compare_dir,"all_libs.csv"),stringsAsFactors = FALSE)


### create interleaved gapfilled file of all markers for tree building
aln_files = list.files(align_dir,pattern=aln_suffix)
print(aln_files)
genelist= read.table(file.path(data_dir,"genelists",filter_file),stringsAsFactors = FALSE)
genelist=genelist$V1
print(genelist)
print(paste("last gene:",genelist[length(genelist)]))

genes = list()

gene_pos_csvfile=file.path(out_dir,paste0(outfile_basename,"_gene_pos.csv"))
gene_pos_nexfile=file.path(out_dir,paste0(outfile_basename,"_gene_pos.nex"))

write("gene,start,end",file=gene_pos_csvfile)
write("#nexus",file=gene_pos_nexfile)
write("begin sets;",file=gene_pos_nexfile,append=TRUE)

charpartition_str="[charpartition all_loc = "  ## [ is to comment out the charpartition string since it doesn't play well with iqtree

start=1
n=1

for(file in aln_files){
  print(file)

  ##check if gene passed filter
  gene=substr(file,1,4)
  if(!(gene %in% genelist) ){
    print("rejected")
    next
  }

  aln = ape::read.dna(file.path(align_dir,file),format=aln_format)

  # check if labels need to be renamed form sample id to species
  if(any(grepl("V",dimnames(aln)[[1]]))){
	print("labels already renamed to species")
  } else{
	dimnames(aln)[[1]] = gsub(paste0("-",gene),"",dimnames(aln)[[1]])

	# create data frame of sample-species pairs in same order as alignment, format species name
 	species = libs[which(libs$Sample %in% dimnames(aln)[[1]]),c("Sample","Species")]
	species =  species[match(dimnames(aln)[[1]],species$Sample),]
	species$Species=gsub(" ","_",species$Species)

	# rename alignment sample IDs as species
  	dimnames(aln)[[1]] = species$Species
 	print(dimnames(aln)[[1]])
  }

  #add to list
  genes[[file]] = aln

  ## print gene length to charset list for nexus file
  # add length of gene to startpoint
  end=start+dim(aln)[[2]]-1

  write(paste(gene,start,end,sep=","),file=gene_pos_csvfile,append=TRUE)
  write(paste0("    charset gene",gene,"  = ",start,"-",end,";"),file=gene_pos_nexfile,append=TRUE)

  ## update charpartition string
  if(gene == genelist[length(genelist)]){
    ## last gene, add ; at end of string
    charpartition_str=paste0(charpartition_str," gene",gene,":",start,"-",end,";]")
  } else {
    charpartition_str=paste0(charpartition_str," gene",gene,":",start,"-",end,",")
  }

  # update start for next loop
  n=n+1
  start=end+1
}
write(paste0("\n    ",charpartition_str,"\nend;"),file=gene_pos_nexfile,append=TRUE)

print("Concatenating genes")
all_markers = genes[[1]]
for(i in 2:length(genes)){
  print(names(genes)[i])
  all_markers = cbind(all_markers,genes[[i]],fill.with.gaps=TRUE)
}

phy_outfn=file.path(out_dir,paste0(outfile_basename,".phy"))
print(paste0("Writing concatenated alignments in phylip format: ",phy_outfn))
write.dna(all_markers,file=phy_outfn,format="interleaved") ##the input RAxML wants

nex_outfn=file.path(out_dir,paste0(outfile_basename,".nex"))
print(paste0("Writing concatenated alignments in nexus format: ",nex_outfn))
write.nexus.data(all_markers,file=nex_outfn,datablock=FALSE)
## for IQtree ModelFinder, comment out the charpartition line with []
##for PAUP (SVDquartets); add `Set allowPunct=yes;` to top
