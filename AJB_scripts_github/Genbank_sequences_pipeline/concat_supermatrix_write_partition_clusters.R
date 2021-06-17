if(!require(ape)) install.packages("ape",repos="https://cran.ma.imperial.ac.uk/")
if(!require(ape)) install.packages("dplyr",repos="https://cran.ma.imperial.ac.uk/")
library(ape)
library(dplyr)

## arguments: alignment directory name, output concat file basename, filter filename (eg supercontigs_secondpass.txt)

data_dir = "~/HybPiper/genbank_compare/"

align_dir=file.path(data_dir,"alignments/trimmed")
outfile_basename="Veronica_genbank_concat5"

print(align_dir)
print(outfile_basename)

aln_suffix="*trm*"
aln_format="fasta"



### create interleaved gapfilled file of all markers for IQTree
aln_files = list.files(align_dir,pattern=aln_suffix)
print(aln_files)

genes = list()

gene_pos_csvfile=file.path(data_dir,paste0(outfile_basename,"_gene_pos.csv"))
gene_pos_nexfile=file.path(data_dir,paste0(outfile_basename,"_gene_pos.nex"))

write("gene,start,end",file=gene_pos_csvfile)
write("begin sets;",file=gene_pos_nexfile)

charpartition_str="charpartition all_loc = "

start=1
n=1

for(file in aln_files){
  print(file)
  
 gene=gsub("Veronica_genbank_","",file)
 gene=gsub("_aln_trm.fasta","",gene)
    
  aln = ape::read.dna(file.path(align_dir,file),format=aln_format)
  ## get rid of accession number in names
  attributes(aln)$dimnames[[1]]= gsub("^[[:alnum:]]*_","",attributes(aln)$dimnames[[1]])

  #add to list
  genes[[file]] = aln

  ## print gene length to charset list for nexus file
  # add length of gene to startpoint   
  end=start+dim(aln)[[2]]-1

  write(paste(gene,start,end,sep=","),file=gene_pos_csvfile,append=TRUE)
  write(paste0("    charset ",gene,"  = ",start,"-",end,";"),file=gene_pos_nexfile,append=TRUE)

  ## update charpartition string
  if(n < length(aln_files)){
    charpartition_str=paste0(charpartition_str," ",gene,":",start,"-",end,",")
  } else {
    charpartition_str=paste0(charpartition_str," ",gene,":",start,"-",end,";")
  }
  
  # update start for next loop
  n=n+1
  start=end+1
}
write(paste0("\n    ",charpartition_str,"\nend;"),file=gene_pos_nexfile,append=TRUE)

print("Concatenating genes")
all_markers = genes[[1]]
for(i in 2:length(genes)){
  all_markers = cbind(all_markers,genes[[i]],fill.with.gaps=TRUE)
}

print(paste0("Writing concatenated alignments in phylip format: ",data_dir,outfile_basename,".phy"))
write.dna(all_markers,file=file.path(data_dir,paste0(outfile_basename,".phy")),format="interleaved") ##the input RAxML wants

print(paste0("Writing concatenated alignments in nexus format: ", data_dir,outfile_basename,".nex"))
write.nexus.data(all_markers,file=file.path(data_dir,paste0(outfile_basename,".nex")),datablock=FALSE) ##for PAUP (SVDquartets); add `Set allowPunct=yes;` to top

