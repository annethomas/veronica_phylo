#######################################################################
### after trimming: filter alignments based on length and %identity ###
#######################################################################
library(dplyr)
library(ape)

hyb_dir="/production/at820/genbank_subset"
compare_dir="~/HybPiper/genbank_compare"



## input_dir should be relative filepath only with no trailing /
input_dir=commandArgs(trailingOnly = TRUE)[1] ## e.g. exons_inframe_gsubset_aln

if(is.na(input_dir)){
	stop(paste("please provide input alignment directory name --relative path, no '/'; assuming it is in",hyb_dir))
}

aln_dir=file.path(hyb_dir,input_dir)
stats_dir=file.path(compare_dir,"stats",paste0(input_dir,"_stats"))
dir.create(stats_dir)

file_basename=input_dir

## identity stats step: get_alignment_identity_stats
#print("running get_alignment_identity_stats.sh")
system(paste("~/HybPiper/genbank_compare/scripts/get_alignment_identity_stats.sh",input_dir))


print("compiling additional alignment stats")
system(paste("~/HybPiper/genbank_compare/scripts/get_gene_stats.sh", input_dir))

## PIS

aln_files = list.files(aln_dir,pattern="*.fasta")
pis_df=data.frame(gene=character(0),pis=numeric(0),pis_pct=numeric(0))

for(file in aln_files){
  gene=substr(file,1,4)
  print(gene)
  print(file.path(aln_dir,file))

  aln = ape::read.dna(file.path(aln_dir,file),format="fasta")
  aln_length = dim(aln)[2]
  print(aln_length)
  pis=phyloch::pis(aln)
  pis_pct=pis/dim(aln)[2]

  pis_df=rbind(pis_df,data.frame(gene=gene,pis=pis,pis_pct=pis_pct))
}




## filtering
print("filtering genes")

ident_thr = .655
occ_thr = 39
len_thr = 150
pis_thr = 20

stats_fn=file.path(stats_dir,paste0(file_basename,"_stats.txt"))
print(stats_fn)
stats=read.table(stats_fn,sep="\t",header=TRUE)
stats$gene=substr(stats$file,1,4)

ident_fn=file.path(stats_dir,paste0(file_basename,"_avg_identity.csv"))
print(ident_fn)
ident=read.csv(ident_fn,stringsAsFactors = FALSE)
ident$gene=as.character(ident$gene)

stats_ident=dplyr::left_join(stats,ident,by="gene")
stats_ident=dplyr::left_join(stats_ident,pis_df,by="gene")

print(head(stats_ident))

## filtering
secondpass=stats_ident[stats_ident$num_seqs > occ_thr &
                            stats_ident$avg_identity > ident_thr &
                            stats_ident$max_len > len_thr &
			    stats_ident$pis >= 20,]

write.csv(stats_ident,file.path(stats_dir,paste0(file_basename,"_secondpass_stats.csv")))

secondpass_out_fn=file.path(hyb_dir,paste0(file_basename,"_secondpass.txt"))
write.table(as.integer(secondpass$gene),secondpass_out_fn,row.names=FALSE,col.names=FALSE)

print(paste(length(secondpass$gene),"secondpass genes printed to", secondpass_out_fn))
