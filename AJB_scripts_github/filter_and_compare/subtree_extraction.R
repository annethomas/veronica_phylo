###############
## subgroups ##
###############
### print group-colored trees
## prep species info
libs_info=read.csv(file.path(hyb_dir,"libs_info_20191124.csv"),stringsAsFactors=FALSE)
libs_info$Species=gsub("V. ","Veronica_",libs_info$Species)
libs_info[which(libs_info$Species=="Veronica_vernicosa-2926"),"Species"]="Veronica_vernicosa"
libs_info[which(libs_info$Species=="Veronica_hookeriana-BDM12"),"Species"]="Veronica_hookeriana"

full_tree_names=c(exons="Veronica_173exons_gsubset_concat",
                  introns="Veronica_167introns_gsubset_concat",
                  supercontigs="Veronica_229supercontigs_gsubset_concat",
                  genbank="Veronica_genbank_concat5")
intersect_tree_names=c(exons="Veronica_exons_intersect_concat",
                       introns="Veronica_introns_intersect_concat",
                       supercontigs="Veronica_supercontigs_intersect_concat")
sortadate_strict_tree_names=c(exons="Veronica_med51exons_gsubset_concat",
                              introns="Veronica_med63introns_gsubset_concat",
                              supercontigs="Veronica_med80supercontigs_gsubset_concat")
sortadate_lax_tree_names=c(exons="Veronica_fix96exons_gsubset_concat",
                           introns="Veronica_fix116introns_gsubset_concat",
                           supercontigs="Veronica_fix150supercontigs_gsubset_concat")  

all_nodekeys=read.csv(file.path(stat_dir,"compare_datasets_nodekey.csv"))


## find nodes of groups
## load add_node_group
source('C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/subtree_functions.R')

#genbank_species=unlist(read.table(file.path(compare_dir,"genbank_compare_species.txt"),stringsAsFactors = FALSE))

# set which dataset to use
nodekeys=full_nodekeys
tree_names=full_tree_names
rooted_tree_filepaths=file.path(compare_dir,"trees",paste0(full_tree_names,"_aperr.tre"))
names(rooted_tree_filepaths)=names(tree_names)

## create tailored list of groups to allow for splitting up paraphyletic groups
species_info=libs_info[,c("Species","hebe_group","hebe_group_detail")]

species_info[which(species_info$hebe_group=="whipcord_hebes"),"hebe_group"]="hebes"

species_info[which(species_info$Species=="Veronica_macrantha"),"hebe_group"]="hebes?"
species_info[which(species_info$Species=="Veronica_cupressoides"),"hebe_group"]="semi-whipcord_hebes?"

species_info[which(species_info$Species %in% c("Veronica_senex","Veronica_lilliputiana")),"hebe_group_detail"]="speedwell_hebes2"
species_info[which(species_info$hebe_group_detail=="speedwell_hebes"),"hebe_group_detail"]="speedwell_hebes1"

species_info[which(species_info$Species %in% c("Veronica_armstrongii","Veronica_salicornioides","Veronica_ochracea")),"hebe_group_detail"]="whipcord_hebes1"
species_info[which(species_info$hebe_group_detail=="whipcord"),"hebe_group_detail"]="whipcord_hebes2"

### supercontigs
cf_stree_rr=ape::read.tree(rooted_tree_filepaths["supercontigs"])
cf_stree_nodekey=nodekeys[["supercontigs"]]


## monophyletic
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"sun_hebes")
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"snow_hebes")
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"semi-whipcord_hebes")
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"speedwell_hebes1")
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"hebes")
# subgroups
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"whipcord_hebes1",subgroup=TRUE)
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"whipcord_hebes2",subgroup=TRUE)
cf_stree_nodekey=add_node_group(cf_stree_rr,cf_stree_nodekey,species_info,"buxifoliatae",subgroup=TRUE)

write.csv(cf_stree_nodekey,file.path(stat_dir,"subgroups_supercontig_nodekey.csv"))
## stats
cf_stree_nodekey_gr=group_by(cf_stree_nodekey,group, subgroup)
cf_stree_nodekey_gr=group_by(cf_stree_nodekey,group)

cf_stree_group_stats=summarize(cf_stree_nodekey_gr,count.nodes=n(),min.bs=min(bootstrap),mean.bs=mean(bootstrap),max.bs=max(bootstrap),
                               min.gcf=min(gCF),mean.gcf=mean(gCF),max.gCF=max(gCF),min.scf=min(sCF),mean.sCF=mean(sCF),max.scf=max(sCF))
ggplot(data=cf_stree_nodekey_gr,aes(x=group,y=bootstrap))+
  geom_boxplot()
ggplot(data=cf_stree_nodekey_gr,aes(x=group,y=gCF))+
  geom_boxplot()
ggplot(data=cf_stree_nodekey_gr,aes(x=group,y=sCF))+
  geom_boxplot()

### exons
cf_etree_rr=ape::read.tree(rooted_tree_filepaths["exons"])
cf_etree_nodekey=nodekeys[["exons"]]

species_info_exons=species_info
species_info_exons[which(species_info_exons$Species == "Veronica_lilliputiana"),"hebe_group_detail"]="speedwell_hebes1"

## monophyletic
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,
                                species_info_exons,"sun_hebes")
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,
                                species_info_exons,"snow_hebes")
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,
                                species_info_exons,"semi-whipcord_hebes")
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,
                                species_info_exons,"speedwell_hebes1")
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,
                                species_info_exons,"hebes")
# subgroups
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"whipcord_hebes1",subgroup=TRUE)

cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"whipcord_hebes2",subgroup=TRUE)
cf_etree_nodekey=add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"buxifoliatae",subgroup=TRUE)

## stats
cf_etree_nodekey_gr=group_by(cf_etree_nodekey,group, subgroup)
cf_etree_nodekey_gr=group_by(cf_etree_nodekey,group)

cf_etree_group_stats=summarize(cf_etree_nodekey_gr,count.nodes=n(),min.bs=min(bootstrap),mean.bs=mean(bootstrap),max.bs=max(bootstrap),
                               min.gcf=min(gCF),mean.gcf=mean(gCF),max.gCF=max(gCF),min.scf=min(sCF),mean.sCF=mean(sCF),max.scf=max(sCF))
ggplot(data=cf_etree_nodekey_gr,aes(x=group,y=bootstrap))+
  geom_boxplot()
ggplot(data=cf_etree_nodekey_gr,aes(x=group,y=gCF))+
  geom_boxplot()
ggplot(data=cf_etree_nodekey_gr,aes(x=group,y=sCF))+
  geom_boxplot()

#paraphyletic
add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"apertae_large",subgroup=TRUE)
add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"subcarnosae",subgroup=TRUE)
add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"apertae_small",subgroup=TRUE)
add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"occlusae",subgroup=TRUE)
add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"connatae",subgroup=TRUE)
#add_node_group(cf_etree_rr,cf_etree_nodekey,species_info_exons,"grandiflorae") #only one species

## introns
cf_itree_rr=ape::read.tree(rooted_tree_filepaths["introns"])
cf_itree_nodekey=nodekeys[["introns"]]

species_info_introns=species_info
#species_info_introns[which(species_info_introns$hebe_group=="speedwell_hebes"),"hebe_group_detail"]="speedwell_hebes1"
species_info_introns[which(species_info_introns$Species == "Veronica_lilliputiana"),
                     "hebe_group_detail"]="speedwell_hebes3"
species_info_introns[which(species_info_introns$Species %in% 
                             c("Veronica_colostylis","Veronica_plano-petiolata","Veronica_linifolia")),
                     "hebe_group_detail"]="speedwell_hebes2"

## monophyletic
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"sun_hebes")
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"snow_hebes")
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"semi-whipcord_hebes")
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"speedwell_hebes1")
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"hebes")
## check:
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,
                                species_info_introns,"speedwell_hebes2")
# subgroups
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"whipcord_hebes1",subgroup=TRUE)

cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"whipcord_hebes2",subgroup=TRUE)
cf_itree_nodekey=add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"buxifoliatae",subgroup=TRUE)

## stats
cf_itree_nodekey_gr=group_by(cf_itree_nodekey,group, subgroup)
cf_itree_nodekey_gr=group_by(cf_itree_nodekey,group)

cf_itree_group_stats=summarize(cf_itree_nodekey_gr,count.nodes=n(),min.bs=min(bootstrap),mean.bs=mean(bootstrap),max.bs=max(bootstrap),
                               min.gcf=min(gCF),mean.gcf=mean(gCF),max.gCF=max(gCF),min.scf=min(sCF),mean.sCF=mean(sCF),max.scf=max(sCF))
ggplot(data=cf_itree_nodekey_gr,aes(x=group,y=bootstrap))+
  geom_boxplot()
ggplot(data=cf_itree_nodekey_gr,aes(x=group,y=gCF))+
  geom_boxplot()
ggplot(data=cf_itree_nodekey_gr,aes(x=group,y=sCF))+
  geom_boxplot()

#paraphyletic
add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"apertae_large",subgroup=TRUE)
add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"subcarnosae",subgroup=TRUE)
add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"apertae_small",subgroup=TRUE)
add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"occlusae",subgroup=TRUE)
add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"connatae",subgroup=TRUE)
#add_node_group(cf_itree_rr,cf_itree_nodekey,species_info_introns,"grandiflorae") #only one species


## genbank
gtree_rr=ape::read.tree(rooted_tree_filepaths["genbank"])
gtree_nodekey=full_nodekeys[["genbank"]]

species_info_genbank=species_info
species_info_genbank[which(species_info_genbank$hebe_group=="speedwell_hebes"),"hebe_group_detail"]="speedwell_hebes1"
species_info_genbank[which(species_info_genbank$Species == "Veronica_plano-petiolata"),
                     "hebe_group_detail"]="speedwell_hebes2"


## monophyletic
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
                             species_info_genbank,"sun_hebes")
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
                             species_info_genbank,"snow_hebes")
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
                             species_info_genbank,"semi-whipcord_hebes")
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
                             species_info_genbank,"speedwell_hebes1")
# gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
#                                 species_info_genbank,"speedwell_hebes2")
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,
                             species_info_genbank,"hebes")
# subgroups (haven't resolved)
gtree_nodekey$subgroup=rep(NA,nrow(gtree_nodekey))
gtree_nodekey$subgroup_mrca=rep(FALSE,nrow(gtree_nodekey))

gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,species_info_genbank,"whipcord_hebes1",subgroup=TRUE)
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,species_info_genbank,"whipcord_hebes2",subgroup=TRUE)
gtree_nodekey=add_node_group(gtree_rr,gtree_nodekey,species_info_genbank,"buxifoliatae",subgroup=TRUE)

## stats
#gtree_nodekey_gr=group_by(gtree_nodekey,group, subgroup)
gtree_nodekey_gr=group_by(gtree_nodekey,group)

gtree_group_stats=summarize(gtree_nodekey_gr,count.nodes=n(),min.bs=min(bootstrap),mean.bs=mean(bootstrap))
ggplot(data=gtree_nodekey_gr,aes(x=group,y=bootstrap))+
  geom_boxplot()

#### combined boxplots ######
all_trees_nodekey=rbind(cf_stree_nodekey,cf_etree_nodekey,cf_itree_nodekey,gtree_nodekey)

## add columns with number of taxa in group
group_counts=all_trees_nodekey %>% group_by(tree_type,group) %>% summarize(count=n())
all_trees_nodekey=left_join(all_trees_nodekey,group_counts,by=c("tree_type","group"))
all_trees_nodekey$mrca_bs=ifelse(all_trees_nodekey$mrca,all_trees_nodekey$bootstrap,NA)

write.csv(all_trees_nodekey,file.path(stat_dir,"compare_subgroups_nodekey_full.csv"),row.names=FALSE)

all_trees_nodekey=read.csv(file.path(stat_dir,"compare_subgroups_nodekey_full.csv"))


all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=group,y=bootstrap))+
  geom_boxplot(aes(fill=tree_type), varwidth=TRUE) #+
#  theme(axis.text.x = element_text(angle = 90))
## specify colors to match the three-subset figs
## add count as text annotation

### group support point/linerange quantiles with size scale
all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=group,y=bootstrap))+
  geom_linerange(mapping = aes(color=tree_type,size=count),position=position_dodge(width=.7),
                 stat = "summary",
                 linetype="dotted",
                 fun.ymin = function(z) {quantile(z,0.05)},
                 fun.ymax = function(z) {quantile(z,0.95)}) +
  geom_pointrange(mapping = aes(color=tree_type,size=count),position=position_dodge(width=.7),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median) +
  scale_size(range = c(.5, 2),breaks=c(2,10,40)) +
  scale_color_manual(values=c("gray30","gray50","gray90","black"))

all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  ggplot(aes(x=group,y=bootstrap))+
  geom_boxplot(aes(fill=tree_type))

all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  ggplot(aes(x=group,y=gCF))+
  geom_boxplot(aes(fill=tree_type))

all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  ggplot(aes(x=group,y=sCF))+
  geom_boxplot(aes(fill=tree_type))

### node age
all_nodekeys=read.csv(file=file.path(stat_dir,"sortadate/compare_datasets_nodekey.csv"), stringsAsFactors = FALSE)

all_nodekeys %>%
  mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=bootstrap))+
  geom_boxplot(aes(fill=tree_type))

all_nodekeys %>% filter(tree_type!="genbank") %>%
  #mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=gCF))+
  geom_boxplot(aes(fill=tree_type))

all_nodekeys %>% filter(tree_type!="genbank") %>%
  #mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=sCF))+
  geom_boxplot(aes(fill=tree_type))


### extracting actual subtrees
library(TreeDist)

get_subgroup=function(tree, tree_type_str, subgroup_str){
  this_mrca=all_trees_nodekey %>% 
    filter(tree_type==tree_type_str & group==subgroup_str & mrca) %>%
    select(node) %>% unlist()
  print(this_mrca)
  gr=extract.clade(tree,this_mrca)
  plot(gr)
  print(length(gr$tip.label))
  return(gr)
}

unique(all_trees_nodekey[which(all_trees_nodekey$tree_type=="supercontigs"),"group"])

cf_stree_rr=ape::read.tree(rooted_tree_filepaths["supercontigs"])
hebes_s=get_subgroup(cf_stree_rr,"supercontigs","hebes")
speedwell_s=get_subgroup(cf_stree_rr,"supercontigs","speedwell_hebes1")
sun_s=get_subgroup(cf_stree_rr,"supercontigs","sun_hebes")
snow_s=get_subgroup(cf_stree_rr,"supercontigs","snow_hebes")
semi_whip_s=get_subgroup(cf_stree_rr,"supercontigs","semi-whipcord_hebes")

cf_itree_rr=ape::read.tree(rooted_tree_filepaths["introns"])
hebes_i=get_subgroup(cf_itree_rr,"introns","hebes")
speedwell_i=get_subgroup(cf_itree_rr,"introns","speedwell_hebes1")
sun_i=get_subgroup(cf_itree_rr,"introns","sun_hebes")
snow_i=get_subgroup(cf_itree_rr,"introns","snow_hebes")
semi_whip_i=get_subgroup(cf_itree_rr,"introns","semi-whipcord_hebes")

cf_etree_rr=ape::read.tree(rooted_tree_filepaths["exons"])
hebes_e=get_subgroup(cf_etree_rr,"exons","hebes")
speedwell_e=get_subgroup(cf_etree_rr,"exons","speedwell_hebes1")
sun_e=get_subgroup(cf_etree_rr,"exons","sun_hebes")
snow_e=get_subgroup(cf_etree_rr,"exons","snow_hebes")
semi_whip_e=get_subgroup(cf_etree_rr,"exons","semi-whipcord_hebes")

cf_gtree_rr=ape::read.tree(rooted_tree_filepaths["genbank"])
hebes_g=get_subgroup(cf_gtree_rr,"genbank","hebes")
speedwell_g=get_subgroup(cf_gtree_rr,"genbank","speedwell_hebes1")
sun_g=get_subgroup(cf_gtree_rr,"genbank","sun_hebes")
snow_g=get_subgroup(cf_gtree_rr,"genbank","snow_hebes")
semi_whip_g=get_subgroup(cf_gtree_rr,"genbank","semi-whipcord_hebes")

TreeDistance(cf_stree_rr,cf_itree_rr)
TreeDistance(list(hebes_g,e=hebes_e,i=hebes_i,s=hebes_s))
#TreeDistance(speedwell_s,speedwell_i)
TreeDistance(list(sun_g,sun_e,sun_i,sun_s))
TreeDistance(list(snow_g,snow_e,snow_i,snow_s))
#TreeDistance(list(semi_whip_e,semi_whip_i,semi_whip_s))
