#######################################################################
### after trimming: filter alignments based on length and %identity ###
#######################################################################

##Rscript ~/HybPiper/genbank_compare/scripts/create_secondpass_genelist.R exons_gsubset ~/HybPiper/genbank_compare/stats/trm_stats/exons_inframe_gsubset_trm_stats.txt ~/HybPiper/genbank_compare/stats/trm_stats/exons_inframe_gsubset_trm_avg_identity.csv

library(dplyr)

hyb_dir="/production/at820"
compare_dir="~/HybPiper/genbank_compare"

stats_dir=file.path(hyb_dir,"stats")

seq_type=commandArgs(trailingOnly = TRUE)[1] ## supercontig, exon, exons_gsubset, etc
stats_fn=commandArgs(trailingOnly = TRUE)[2] ## ~/HybPiper/genbank_compare/stats/trm_stats/exons_inframe_gsubset_trm_stats.txt
ident_fn=commandArgs(trailingOnly = TRUE)[3] ## ~/HybPiper/genbank_compare/stats/trm_stats/exons_inframe_gsubset_trm_avg_identity.csv

secondpass_out_fn = file.path(compare_dir,paste0(seq_type,"_secondpass.txt"))

ident_thr = .655
occ_thr = 39
len_thr = 150

ident=read.csv(ident_fn,stringsAsFactors = FALSE)
ident$gene=as.character(ident$gene)

stats=read.table(stats_fn,sep="\t",header=TRUE)
stats$gene=substr(stats$file,1,4)

stats_ident=dplyr::left_join(stats,ident,by="gene")
print(head(stats_ident))

## filtering
secondpass=stats_ident[stats_ident$num_seqs > occ_thr & 
                            stats_ident$avg_identity > ident_thr & 
                            stats_ident$max_len > len_thr,]
write.csv(stats_ident,file.path(compare_dir,paste0(seq_type,"_secondpass_stats.csv")))

write.table(as.integer(secondpass$gene),secondpass_out_fn,row.names=FALSE,col.names=FALSE)





##Possibility: keep exons between .6 and .655 identity if greater than 500 bp or some other length threshold




