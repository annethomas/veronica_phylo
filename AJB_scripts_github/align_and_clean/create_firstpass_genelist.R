#####################################################################################
### After running HybPiper: filter raw genes based on occupancy and paralog detection ###
#####################################################################################
## NB: run interactively

hyb_dir="~/HybPiper/genbank_compare"
stats_dir=file.path(hyb_dir,"stats")

######################################################################################
### original lengths_stats: read in stats and seqlengths from HybPiper and collate ###
######################################################################################
v_stats=read.table(file.path(stats_dir,"hybpiper_stats.txt"),header=TRUE)
v_lengths = read.table(file.path(stats_dir,"gene_lengths.txt"),header=TRUE,row.names = 1)

v_lengths_stats=v_lengths

## add gene summary stats to seqlengths
last_samp=nrow(v_lengths_stats)
v_lengths_stats["gene_length_med",] = apply(v_lengths_stats[2:last_samp,],2,median)
v_lengths_stats["percent_0s",] = apply(v_lengths_stats[2:last_samp,],2,function(x) length(which(x==0))/(last_samp-1))
v_lengths_stats=round(v_lengths_stats[],2)
length(which(v_lengths_stats["percent_0s",]>.5))
zeros_50 = names(v_lengths_stats)[which(v_lengths_stats["percent_0s",]>.5)]

## extra stats
v_lengths_stats["gene_length_mean",] = colMeans(v_lengths_stats[2:last_samp,])
v_lengths_stats["gene_length_cv",] = apply(v_lengths_stats[2:last_samp,],2,cv)
plot(as.numeric(v_lengths_stats["col_cv",]),as.numeric(v_lengths_stats["col_med",]))


## save collated version
write.csv(v_lengths_stats,file=file.path(stats_dir,"lengths_stats.csv"))

######################################################################

###############################################
### read in lengths_stats when running again ##
###############################################
v_lengths_stats = read.csv(lengths_stats_fn,stringsAsFactors = FALSE)
zeros_50 = names(v_lengths_stats)[which(v_lengths_stats["percent_0s",]>.5)]
zeros_50=gsub("X","",zeros_50)

##############################################################################
### original paralog list from paralog_investigate_loop.sh -- extract count ##
##############################################################################

paralog_warnings=read.table(file.path(stats_dir,"paralog_warnings.txt"))
has_paralogs=unique(paralog_warnings$V5)
has_paralogs_count=as.data.frame(table(paralog_warnings$V5))
hist(has_paralogs_count$Freq)
length(which(has_paralogs_count$Freq>11))
write.csv(has_paralogs_count,file.path(stats_dir,"has_paralogs_count.csv"))


###############################################
### read in paralog count when running again ##
###############################################

has_paralogs_count = read.csv(paralogs_count_fn,stringsAsFactors = FALSE)
## genes with paralogs detected in more than 10% of samples (>11/117)
high_paralog = as.character(has_paralogs_count[which(has_paralogs_count$Freq>11),"Var1"])

#############
### filter ##
#############

genes_filter_firstpass=union(zeros_50,high_paralog)

genes=gsub("X","",names(v_lengths_stats))
genes_firstpass=as.numeric(sort(setdiff(genes,genes_filter_firstpass)))

write.table(genes_firstpass,file.path(stats_dir,"veronica_genes_firstpass.txt"),row.names = FALSE)
