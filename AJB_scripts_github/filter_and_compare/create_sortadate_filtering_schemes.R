library(ggplot2)
library(viridis)
library(dplyr)
library(ape)
library(phytools)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
hyb_dir=file.path(data_dir,"HybPiper")
compare_dir=file.path(data_dir,"genbank_compare")
stat_dir=file.path(compare_dir,"stats")

###############
## sortadate ##
###############


## everything rooted
exons_sorted=read.table(file.path(stat_dir,"sortadate/exons_gsubset_trm_cleaned_tiptrim_sortedgenes.rooted.txt"),header=TRUE)
exons_sorted$sort=1:nrow(exons_sorted)
exons_sorted$gene=substr(exons_sorted$name,1,4)

compare=dplyr::left_join(exons_sortedu,exons_sorted,by="name")
boxplot(exons_sortedu$bipartition,exons_sorted$bipartition)

sups_sorted=read.table(file.path(stat_dir,"sortadate/supercontigs_gsubset_trm_cleaned_tiptrim_sortedgenes.txt"),header=TRUE)
sups_sorted$sort=1:nrow(sups_sorted)
sups_sorted$gene=substr(sups_sorted$name,1,4)

introns_sorted=read.table(file.path(stat_dir,"sortadate/introns_gsubset_trm_cleaned_tiptrim_sortedgenes.txt"),header=TRUE)
introns_sorted$sort=1:nrow(introns_sorted)
introns_sorted$gene=substr(introns_sorted$name,1,4)


boxplot(exons_sorted$bipartition,introns_sorted$bipartition,sups_sorted$bipartition)
boxplot(exons_sorted$treelength,introns_sorted$treelength,sups_sorted$treelength)
boxplot(exons_sorted$root.to.tip_var,introns_sorted$root.to.tip_var,sups_sorted$root.to.tip_var)

summary(sups_sorted$root.to.tip_var)

filter_genes=function(sorted_list,bipart_thr,length_thr,var_thr,quantile=FALSE){
  
  if(!quantile){
    return(dplyr::filter(sorted_list,bipartition>bipart_thr,
                         treelength>length_thr,
                         root.to.tip_var<var_thr))
  }
  else{
    return(dplyr::filter(sorted_list,bipartition>=quantile(sorted_list$bipartition,bipart_thr),
                         treelength>=quantile(sorted_list$treelength,length_thr),
                         root.to.tip_var<=quantile(sorted_list$root.to.tip_var,var_thr)))
  }
}


exon_quant1=filter_genes(exons_sorted,.25,.25,1,quantile=TRUE)
#nrow(filter_genes(exons_sorted,.5,0,1,quantile=TRUE))
exon_quant2=filter_genes(exons_sorted,.5,0,1,quantile=TRUE)

exon_quant=filter_genes(exons_sorted,.5,.5,1,quantile=TRUE)
median(exons_sorted$bipartition) #0.02597403
median(exons_sorted$treelength) # 0.6385135
min(exons_sorted$treelength) #0.224922
nrow(exon_quant) #51
exon_fixed=filter_genes(exons_sorted,.02,.1,.024,quantile=FALSE)
nrow(exon_fixed) #96 (same bp as quantile but with no low-length removed)

intron_quant=filter_genes(introns_sorted,.5,.5,1,quantile=TRUE)
nrow(intron_quant) #63
median(introns_sorted$bipartition) #0.02597403
median(introns_sorted$treelength) #1.10881
min(introns_sorted$treelength) 0.411229
intron_fixed=filter_genes(introns_sorted,.02,.5,.024,quantile=FALSE)
nrow(intron_fixed) #116 (same bp as quantile but with fewer low-length removed)

sups_quant=filter_genes(sups_sorted,.5,.5,1,quantile=TRUE)
median(sups_sorted$bipartition) #0.03896104
median(sups_sorted$treelength) # 0.792527
min(sups_sorted$treelength) #0.22995
nrow(sups_quant) #80
sups_fixed=filter_genes(sups_sorted,.02,.5,.024,quantile=FALSE)
nrow(sups_fixed) #150 (more of both low-bp and low-length)

write.csv(exon_quant,file=file.path(stat_dir,"exons_med51.csv"))
write.csv(exon_fixed,file=file.path(stat_dir,"exons_fix96.csv"))
write.table(sort(as.numeric(exon_quant$gene)),
            file=file.path(stat_dir,"genelists/exons_med51_genelist.txt"),
            row.names=FALSE,col.names=FALSE)
write.table(sort(as.numeric(exon_fixed$gene)),
            file=file.path(stat_dir,"genelists/exons_fix96_genelist.txt"),
            row.names=FALSE,col.names=FALSE)


write.csv(intron_quant,file=file.path(stat_dir,"introns_med63.csv"))
write.csv(intron_fixed,file=file.path(stat_dir,"introns_fix116.csv"))
write.table(sort(as.numeric(intron_quant$gene)),
            file=file.path(stat_dir,"genelists/introns_med63_genelist.txt"),
            row.names=FALSE,col.names=FALSE)
write.table(sort(as.numeric(intron_fixed$gene)),
            file=file.path(stat_dir,"genelists/introns_fix116_genelist.txt"),
            row.names=FALSE,col.names=FALSE)

write.csv(sups_quant,file=file.path(stat_dir,"supercontigs_med80.csv"))
write.csv(sups_fixed,file=file.path(stat_dir,"supercontigs_fix150.csv"))
write.table(sort(as.numeric(sups_quant$gene)),
            file=file.path(stat_dir,"genelists/supercontigs_med80_genelist.txt"),
            row.names=FALSE,col.names=FALSE)
write.table(sort(as.numeric(sups_fixed$gene)),
            file=file.path(stat_dir,"genelists/supercontigs_fix150_genelist.txt"),
            row.names=FALSE,col.names=FALSE)