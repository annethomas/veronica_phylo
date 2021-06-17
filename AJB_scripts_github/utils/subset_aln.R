if(!require(ape)) install.packages("ape",repos="https://cran.ma.imperial.ac.uk/")
if(!require(ape)) install.packages("dplyr",repos="https://cran.ma.imperial.ac.uk/")
library(ape)
library(dplyr)

## arguments: alignment directory name, output concat file basename, filter filename (eg supercontigs_secondpass.txt)

hyb_dir= "~/HybPiper/"

input_dir=commandArgs(trailingOnly=TRUE)[1] ## where the alignments are
out_dir=commandArgs(trailingOnly=TRUE)[2]
species_filter_file=commandArgs(trailingOnly=TRUE)[3] # full path--species in "sample_V._species" format
gene_filter_file=commandArgs(trailingOnly=TRUE)[4] ## full path; only necessary if only subsetting some of the alignments in input_dir



dir.create(out_dir)

print(out_dir)

aln_suffix=".fasta"
aln_format="fasta"

libs=read.csv(file.path(hyb_dir,"all_libs.csv"),stringsAsFactors = FALSE)
print("test")

specieslist= read.table(species_filter_file,stringsAsFactors = FALSE)
specieslist=specieslist$V1

#specieslist= c("S2L_V._scrupea","S4C_V._pentasepala","H2E_V._lavaudiana","E3_V._hulkeana","H2F_V._maccaskillii","S3D_V._raoulii","A3_V._densifolia")


print(specieslist)

genefilter=FALSE
if(!is.na(gene_filter_file)){
  genefilter=TRUE
  genelist= read.table(gene_filter_file,stringsAsFactors = FALSE)
  genelist=genelist$V1
}


### 
aln_files = list.files(input_dir,pattern=aln_suffix)
print(aln_files)

for(file in aln_files){
  print(file)

  gene=substr(file,1,4)

  if(genefilter){
    ##check if gene passed filter
    if(!(gene %in% genelist) ){
      print("rejected")
      next
    }
  }
  
  aln = ape::read.dna(file.path(input_dir,file),format=aln_format)

  # check if labels need to be renamed from sample id to species
  if(any(grepl("V",dimnames(aln)[[1]]))){
	print("labels already renamed to species")
  } else{
	dimnames(aln)[[1]] = gsub(paste0("-",gene),"",dimnames(aln)[[1]])

	# create data frame of sample-species pairs in same order as alignment, format species name
 	species = libs[which(libs$Sample %in% dimnames(aln)[[1]]),c("Sample","Species")]
	species =  species[match(dimnames(aln)[[1]],species$Sample),]
	species$Species=gsub(" ","_",species$Species)
	species$tip=paste0(species$Sample,"_",species$Species)
	
	# rename alignment sample IDs as species
  	dimnames(aln)[[1]] = species$tip
 	print(dimnames(aln)[[1]])
  }

  aln=aln[which(dimnames(aln)[[1]] %in% specieslist),]
  
  print(dimnames(aln))
  write.dna(aln,file.path(out_dir,file),format="fasta")
  
}

