#if(!require(ape)) install.packages("ape",repos="https://cran.ma.imperial.ac.uk/")
#if(!require(ape)) install.packages("dplyr",repos="https://cran.ma.imperial.ac.uk/")
library(ape)
library(dplyr)

data_dir = "/production/at820/genbank_subset"
input_dir=commandArgs(trailingOnly=TRUE)[1] # e.g. exons_inframe_gsubset_trm
align_dir = file.path(data_dir,input_dir)

filter_fn=commandArgs(trailingOnly=TRUE)[2] #full filepath

aln_suffix="*trm*"
clean_dir=file.path(data_dir,paste0(input_dir,"_cleaned"))
dir.create(clean_dir)

#clean_phy_dir=paste0(data_dir,"exons_trimmed_cleaned_phy/")
#dir.create(clean_phy_dir)


aln_files = list.files(align_dir,pattern=aln_suffix)
filter=unlist(read.table(filter_fn))

for(file in aln_files){
  print(file)
  ##check if gene passed avg identity filter
  gene=substr(file,1,4)
  if(!gene %in% filter){
	next
  }


  aln = ape::read.dna(file.path(align_dir,file),format="fasta")
  orig_len=dim(aln)[[1]]
  
  # delete sequences with less than 5% of the total alignment length
  aln_nogaps=del.rowgapsonly(aln,threshold=.95)
  if(dim(aln_nogaps)[[1]]<orig_len) print(paste("Num rows deleted: ",orig_len-dim(aln_nogaps)[[1]]))

  write.dna(aln_nogaps,file=file.path(clean_dir,file),format="fasta")
    
  #phy_filename=paste0(substr(file,1,nchar(file)-nchar("fasta")),"phy")
  #write.dna(aln_nogaps,file=paste0(clean_phy_dir,phy_filename),format="interleaved") ##the input iqtree wants
}

