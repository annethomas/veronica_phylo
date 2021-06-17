# AJB manuscript plots

library(ggplot2)
library(viridis)
library(dplyr)
library(ape)
library(phytools)
library(ggpubr)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
hyb_dir=file.path(data_dir,"HybPiper")
compare_dir=file.path(data_dir,"genbank_compare")
stat_dir=file.path(compare_dir,"stats")
plot_dir=file.path(stat_dir,"manuscript_plots")
#dir.create(plot_dir)

#################################################################
### final plots and stats ####
#########################################################



#load make_nodekey and compile_nodekeys functions
source('C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/phylogenetics/nodekey_functions.R')

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


#################################################
## concordance factors + bootstrap comparisons ##
#################################################

## print dataset-specific plots and get lists of nodekeys
### accesses rooted (or roots) concordance factor+bootstrap-labelled trees, <tree_name>_aperr.tre in tree_dir
full_nodekeys=compile_nodekeys(full_tree_names, "full",root=FALSE,plots=FALSE)
intersection_nodekeys=compile_nodekeys(intersect_tree_names, "intersection",root=FALSE,plots = FALSE)
sortadate_strict_nodekeys=compile_nodekeys(sortadate_strict_tree_names, "SortaDate_TL",root=FALSE,plots=FALSE)
sortadate_lax_nodekeys=compile_nodekeys(sortadate_lax_tree_names, "SortaDate_BP",root=FALSE,plots=FALSE)

## combine all datasets in one boxplot per metric
full_nodekeys_df=do.call("rbind",full_nodekeys)
intersection_nodekeys=do.call("rbind",intersection_nodekeys)
sortadate_strict_nodekeys_df=do.call("rbind",sortadate_strict_nodekeys)
sortadate_lax_nodekeys_df=do.call("rbind",sortadate_lax_nodekeys)

all_nodekeys=rbind(full_nodekeys_df,intersection_nodekeys,sortadate_strict_nodekeys_df,sortadate_lax_nodekeys_df)
write.csv(all_nodekeys,file=file.path(stat_dir,"compare_datasets_nodekey_updated.csv"),row.names = FALSE)
all_nodekeys=read.csv(file.path(stat_dir,"compare_datasets_nodekey_updated.csv"))

##############
## Figures ###
##############


### FIG 1: photos put together in Inkscape
### FIG 2: methods flowchart in draw.io

### FIG 3: bootstrap/gcf/scf boxplots
#####################################
## for ggarrange panel 
bs=all_nodekeys %>%   
  mutate(tree_type = factor(tree_type, 
                            levels=c("genbank","exons","introns","supercontigs"),
                            labels=c("GenBank","exons","introns","supercontigs"))) %>%  ggplot(aes(x=tree_type,y=bootstrap)) +
  geom_boxplot(aes(fill=dataset),position = position_dodge2(preserve = "single"),size=.25) +
  labs(x="Gene subset",y="Bootstrap (%)",fill="Filtering scheme") +
  theme_classic(base_size=12) +
  scale_fill_brewer(palette="Blues",direction=-1)

gcf=all_nodekeys  %>% filter(tree_type!="genbank") %>%
  ggplot(aes(x=tree_type,y=gCF)) +
  geom_boxplot(aes(fill=dataset),size=.25,outlier.size=0.5) +
  labs(x=NULL,y="gCF (%)",fill="Filtering scheme") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  scale_fill_brewer(palette="Blues",direction=-1)

scf=all_nodekeys %>% filter(tree_type!="genbank") %>%
  ggplot(aes(x=tree_type,y=sCF)) +
  geom_boxplot(aes(fill=dataset),size=.25,outlier.size=0.5)  +
  labs(x="Gene subset",y="sCF (%)",fill="Filtering scheme") +
  theme_classic(base_size=12) +
  scale_fill_brewer(palette="Blues",direction=-1)


gcf_for_legend=all_nodekeys  %>% filter(tree_type!="genbank") %>%
  ggplot(aes(x=tree_type,y=gCF)) +
  geom_boxplot(aes(fill=dataset),size=.25) +
  labs(x=NULL,y="gCF (%)",fill="Filtering scheme") +
  theme_classic(base_size=15) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  scale_fill_brewer(palette="Blues",direction=-1)
legendonly=as_ggplot(get_legend(gcf_for_legend))


### ggarrange
#png(filename = file.path(plot_dir,"filter_compare_box_grid_updated.png"),width=1000,height=750)
pdf(file = file.path(plot_dir,"filter_compare_box_grid_testb.pdf"),width=7.25,height=5.5)

ggpubr::ggarrange(bs,gcf,legendonly,scf,labels=c("(a)","(b)","","(c)"),label.x=-0.02, label.y=1,
                            font.label=list(size=12),legend="none") 
dev.off()

#identical alternative for saving
ggsave(filename="filter_compare_box_grid_testggsave.pdf",path=plot_dir,plot=grid_fig3,device="pdf",
       width=7.25,height=5.5,units="in",dpi=600) # dpi does nothing (this is a vector format) 

#########
### FIG 4: bootstrap, gcf, scf scatterplot
##############################################
### scatterplot of gCF vs sCF with bootstrap and nodetype (full dataset only)

sqrt_it <- scales::trans_new(
  name = "sqrt_scale",
  transform = function (x) { sign(x) * abs(x)^(1/2) },
  inverse = function (x) { sign(x) * abs(x)^2 }
)
pdf(file.path(plot_dir,"gcf_scf_bootstrap_nodeage_fullsuper_scatter-test5.pdf"),width=3.5,height=3.25)

full_nodekeys_df %>% 
  ggplot(aes(x = gCF, y = sCF)) +
  geom_point(aes(shape = parent_of_term, fill = bootstrap), alpha = 0.7, size =1.75) +
  scale_fill_viridis(direction = -1) +
  scale_x_continuous(trans = sqrt_it, limits = c(0, 50), breaks = seq(0, 50, 10)) +
  scale_y_continuous(expand = c(0,1), limits = c(25, 60), breaks = seq(25, 60, 5)) +
  scale_shape_manual(values = c(24, 21), labels  = c("deep", "shallow")) +
  labs(x = "gCF (%)", y = "sCF (%)", shape = "Node type",fill="BS") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",legend.key.height=unit(.1,'in'),legend.key.width=unit(.15,'in'),
        legend.text=element_text(size=8),legend.title=element_text(size=10),
        legend.spacing.x=unit(.025,'in'))
#  theme(legend.position = c(1,0.1), 
#        legend.justification = c(1,0)) 
dev.off()

test=cor.test(full_nodekeys_df$gCF,full_nodekeys_df$bootstrap) 
cor.test.p=function(x,y){
  test=cor.test(x,y)
  return(test$p.value)
}
bs.cf.cors=full_nodekeys_df %>% filter(tree_type!="genbank") %>% group_by(tree_type) %>% 
  summarize(cor.gcf.bs=cor(gCF,bootstrap,use="complete.obs"),
            cor.gcf.bs.p=cor.test.p(gCF,bootstrap),
            cor.scf.bs=cor(sCF,bootstrap,use="complete.obs"),
            cor.scf.bs.p=cor.test.p(sCF,bootstrap))


#####
## FIG 5: supercontig IQtree phylogeny (hand altered in Inkscape)
## plus supplementary tree figures

###############
## full tree ##
###############
### print group-colored trees
## prep species info
#libs_info=read.csv(file.path(hyb_dir,"libs_info_detailed_whipcord.csv"),stringsAsFactors=FALSE)
libs_info=read.csv(file.path(hyb_dir,"libs_info_detailed.csv"),stringsAsFactors=FALSE)
#libs_info$Species=gsub("V. ","Veronica_",libs_info$Species)
libs_info[which(libs_info$Species=="V. vernicosa-2926"),"Species"]="V. vernicosa"
libs_info[which(libs_info$Species=="V. hookeriana-BDM12"),"Species"]="V. hookeriana"


### color by simple hebe group
#hebe_groups = read.csv(file.path(hyb_dir,"hebe_groups_simple_colors.csv"),stringsAsFactors = FALSE)
hebe_groups = read.csv(file.path(hyb_dir,"hebe_groups_simple_colors_gsubset.csv"),stringsAsFactors = FALSE)
hebe_colors = unique(hebe_groups$color)

## IQTREE
#tree_dir=file.path(compare_dir,"trees")
tree_dir=file.path(plot_dir,"trees")
out_dir=tree_dir

make_tree_fig=function(tree_filename,type,legend.pos,scale.pos){
  print(file.path(tree_dir,tree_filename))
  tr=ape::read.tree(file=file.path(tree_dir,tree_filename))
  print("test")
  tr=drop.tip(tr,"Veronica_chamaedrys")
  tr$tip.label=gsub("Veronica_", "V. ",tr$tip.label)
  tips = dplyr::left_join(data.frame("tip"=tr$tip.label),libs_info[,c("Species","Sample","natural_group","hebe_group","hebe_group_detail","alpine_status")],by=c("tip"="Species"))
  tips=dplyr::left_join(tips,hebe_groups)
  
  tr$tip.label=gsub("plano-petiolata", "planopetiolata",tr$tip.label)
  
  ## color nodes above certain bs
  
  
  if(type=="iqtree"){
    
    
    node_labs=tidyr::separate(data.frame(tr$node.label),col=tr.node.label,into=c("bootstrap","gCF","sCF"),sep="/",remove=TRUE,convert=TRUE)
    node_labs$bs_col=ifelse(node_labs$bootstrap>=80,ifelse(node_labs$bootstrap>=99,"turquoise1","grey47"),NA)
    node_labs$node_no=1:nrow(node_labs)+length(tr$tip.label)
    nodes_to_col=node_labs %>% filter(!is.na(bs_col)) %>% dplyr::select(node_no) %>% unlist()
    node_cols=node_labs %>% filter(!is.na(bs_col)) %>% dplyr::select(bs_col) %>% unlist()
    
    legend_specs = read.csv(file.path(plot_dir,"combined_legend_ms.csv"),stringsAsFactors = FALSE)
    
    
    pdf_output_filename = file.path(out_dir,paste0(tree_filename,"_test.pdf"))
    #svg_output_filename = file.path(plot_dir,paste0(tree_filename,"_test.svg"))
    
    pdf(file = pdf_output_filename,width=11,height=15.5)
    plot(tr,tip.color=hebe_colors[tips$code],cex=1,font=2,label.offset = .001)
    nodelabels(tr$node.label,adj = c(-0.1), frame = "n", cex = .8)
    ape::nodelabels(node=nodes_to_col, pch = 19, col = node_cols, cex=1)
    
    legend(x=legend.pos,legend=legend_specs$text,col=legend_specs$color,
           text.font=legend_specs$font,text.col=legend_specs$text.color,
           pch=legend_specs$pch,cex=.8)
    #ape::add.scale.bar(x=0.03,y=0,lwd=2,cex=2)
    ape::add.scale.bar(x=scale.pos[1],y=0,lwd=2,cex=1)
    dev.off()
    
    png_output_filename = file.path(out_dir,paste0(tree_filename,"_test.png"))
    png(filename = png_output_filename,width=2900,height=3900)
    plot(tr,tip.color=hebe_colors[tips$code],cex=3.5,font=2,label.offset = .001)
    nodelabels(tr$node.label,adj = -0.1, frame = "n", cex = 2.5)
    ape::nodelabels(node=nodes_to_col, pch = 19, col = node_cols, cex=4)
    
    legend(x=legend.pos,legend=legend_specs$text,col=legend_specs$color,
           text.font=legend_specs$font,text.col=legend_specs$text.color,
           pch=legend_specs$pch,cex=3.75)
    #ape::add.scale.bar(x=0.03,y=0,lwd=2,cex=2)
    ape::add.scale.bar(x=scale.pos[1],y=0,lwd=2,cex=3)
    dev.off()
  }else {
    
    if(type=="astral"){
      node_labs=tidyr::separate(data.frame(tr$node.label),col=tr.node.label,into=c("LPP","gCF"),sep="/",remove=TRUE,convert=TRUE)
      node_labs$lpp_col=ifelse(node_labs$LPP>.8,ifelse(node_labs$LPP>=.99,"turquoise1","grey47"),NA)
      node_labs$node_no=1:nrow(node_labs)+length(tr$tip.label)
      nodes_to_col=node_labs %>% filter(!is.na(lpp_col)) %>% select(node_no) %>% unlist()
      node_cols=node_labs %>% filter(!is.na(lpp_col)) %>% select(lpp_col) %>% unlist()
      
      
      legend_specs = read.csv(file.path(plot_dir,"combined_legend_astral_ms.csv"),stringsAsFactors = FALSE)
      
    }else{
      node_labs=data.frame(bs=as.numeric(tr$node.label))
      node_labs$bs_col=ifelse(node_labs$bs>80,ifelse(node_labs$bs>=99,"turquoise","gray"),NA)
      node_labs$node_no=1:nrow(node_labs)+length(tr$tip.label)
      nodes_to_col=node_labs %>% filter(!is.na(bs_col)) %>% select(node_no) %>% unlist()
      node_cols=node_labs %>% filter(!is.na(bs_col)) %>% select(bs_col) %>% unlist()
      
      legend_specs = read.csv(file.path(plot_dir,"combined_legend_ms.csv"),stringsAsFactors = FALSE)
    }
    
    if(type=="svdquartets"){
      tr$edge.length = NULL
    }
    
    pdf_output_filename = file.path(out_dir,paste0(tree_filename,"_test.pdf"))
    #svg_output_filename = file.path(plot_dir,paste0(tree_filename,"_test.svg"))
    
    pdf(file = pdf_output_filename,width=11,height=15.5)
    plot(ladderize(tr),tip.color=hebe_colors[tips$code],cex=1,font=2,label.offset = .001)
    nodelabels(tr$node.label,adj = c(1.15,-0.18), frame = "n", cex = .8)
    ape::nodelabels(node=nodes_to_col, pch = 19, col = node_cols, cex=1)
    
    legend(x=legend.pos,legend=legend_specs$text,col=legend_specs$color,
           text.font=legend_specs$font,text.col=legend_specs$text.color,
           pch=legend_specs$pch,cex=.8)
    if(type=="genbank"){
      ape::add.scale.bar(x=scale.pos[1],y=0,lwd=2,cex=1)
    }
    dev.off()
    
    
    png_output_filename = file.path(out_dir,paste0(tree_filename,"_test.png"))
    png(filename = png_output_filename,width=2900,height=3900)
    plot(tr,tip.color=hebe_colors[tips$code],cex=3.5,font=2,label.offset = .001)
    nodelabels(tr$node.label,adj = c(1.15,-0.18), frame = "n", cex = 2.5)
    ape::nodelabels(node=nodes_to_col, pch = 19, col = node_cols, cex=4)
    
    legend(x=legend.pos,legend=legend_specs$text,col=legend_specs$color,
           text.font=legend_specs$font,text.col=legend_specs$text.color,
           pch=legend_specs$pch,cex=3.75)
    if(type=="genbank"){
      ape::add.scale.bar(x=scale.pos[1],y=0,lwd=2,cex=3)
    }
    dev.off()
    
    
  }
  
  
}

make_tree_fig("Veronica_173exons_aperr.tre","iqtree","bottomright",0)
make_tree_fig("Veronica_167introns_aperr.tre","iqtree","bottomright",0)
make_tree_fig("Veronica_229supercontigs_aperr.tre","iqtree","bottomleft",0.055)
make_tree_fig("Veronica_genbank_concat5_aperr.tre","genbank","bottomright",0)

make_tree_fig("Veronica_229supercontigs_gsubset_concat_aperr.astral.tre","astral","bottomleft",0)

make_tree_fig("Veronica_229supercontigs_gsubset.svd.nodes_aperr.tre","svdquartets","bottomleft",0)

###
### FIG 6: line plot of bootstraps by subgroup; only works in  dplyr v 0.8.0.1, ggplot2 v 2.2.1 and below (R v 3.4.1) (see rescue_old_subgroup_plot.R)
### note: resized pdf from old version/computer by printing to pdf at 66% to make it 7.25 inches, then removed margins in Inkscape, but it slightly resizes the lines so they've been shifted slightly
all_trees_nodekey=read.csv(file.path(stat_dir,"compare_subgroups_nodekey_full.csv"),stringsAsFactors=FALSE)

mrca_values=group_by(all_trees_nodekey,dataset,tree_type,group)%>% filter(mrca)

pdf(file=file.path(plot_dir,"subclades_bootstrap_full_pointlinex_gridclassictest.pdf"),width=10,height=8)

all_trees_nodekey %>% filter(!is.na(group),group!="speedwell_hebes2") %>%
  mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  mutate(group = factor(group, levels=c("hebes","speedwell_hebes1","snow_hebes","sun_hebes","semi-whipcord_hebes"),
                        labels=c("hebes (48)","speedwell (10)","snow (6)","sun (5)","semi-whipcord (3)"))) %>%
  ggplot(aes(x=tree_type,y=bootstrap))+
  facet_grid(.~group,switch="x")+
  geom_linerange(mapping = aes(color=tree_type,size=count),position=position_dodge(width=.5),
                 stat = "summary",
                 linetype="dotted",
                 fun.ymin = function(z) {quantile(z,0.05)},
                 fun.ymax = function(z) {quantile(z,0.95)}) +
  geom_pointrange(mapping = aes(color=tree_type,size=count),position=position_dodge(width=.5),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median) +  
  geom_point(aes(x=tree_type,y=mrca_bs,color=tree_type),shape=4,position=position_dodge(width=.5),size=5,stroke=2) +
  #geom_vline(xintercept=c(1.5,2.5,3.5,4.5)) +
  scale_size(range = c(.5, 2),breaks=c(2,10,40)) +
  labs(x="Subclade",y="Bootstrap (%)",color="Gene subset",size="Num. nodes") +
  theme_classic(base_size=16) +
  theme(axis.text.x=element_blank(),axis.ticks = element_blank())+
  theme(legend.position = c(1,0),legend.justification = c(1,0), legend.spacing=unit(0, "inch"),
        #legend.box.background = element_rect(color="black"),
        legend.box.margin = margin(b=6,r=6,t=1))+
  scale_x_discrete(expand=c(0.5,0.5)) +
  scale_color_manual(values=c("gray70","gray50","gray30","black"))
dev.off()

#####
## FIG 7: Consensus network, made in SplitsTree4/Inkscape

#####
## FIG 8: node ages and branch lengths bar plots
######

### summarize branch lengths 
############
stree=read.tree(file.path(compare_dir,"trees/Veronica_229supercontigs_gsubset_concat_aperr.tre"))               
etree=read.tree(file.path(compare_dir,"trees/Veronica_173exons_gsubset_concat_aperr.tre"))               
itree=read.tree(file.path(compare_dir,"trees/Veronica_167introns_gsubset_concat_aperr.tre"))               
gtree=read.tree(file.path(compare_dir,"trees/Veronica_genbank_concat5_aperr.tre"))               


## separate internal and terminal branch lengths
# terminal branches based on the end node in edge[,2]
print_edge_summary=function(tree,tree_type_str){
  term_edge=tree$edge.length[which(tree$edge[,2] %in% 1:length(tree$tip.label))]
  print("terminal edges")
  #print(boxplot(term_edge))
  print(summary(term_edge))
  
  in_edge=tree$edge.length[which(!tree$edge[,2] %in% 1:length(tree$tip.label))]
  print("internal edges")
  #print(boxplot(in_edge))
  print(summary(in_edge))
  
  print(boxplot(term_edge,in_edge))
  print(length(c(in_edge,term_edge)))
  combined_df=data.frame(length=c(in_edge,term_edge),
                         type=c(rep("internal",length(in_edge)),
                                rep("terminal",length(term_edge))),
                         tree_type=tree_type_str)
  return(combined_df)
}

stree=drop.tip(stree,c("Veronica_chamaedrys","Veronica_perfoliata"))
stree_edgelen=print_edge_summary(stree,"supercontigs")

gtree=drop.tip(gtree,c("Veronica_chamaedrys","Veronica_perfoliata"))
gtree_edgelen=print_edge_summary(gtree,"GenBank")

etree=drop.tip(etree,c("Veronica_chamaedrys","Veronica_perfoliata"))
etree_edgelen=print_edge_summary(etree,"exons")

itree=drop.tip(itree,c("Veronica_chamaedrys","Veronica_perfoliata"))
itree_edgelen=print_edge_summary(itree,"introns")

gtree_edgelen$tree_type="GenBank"
etree_edgelen$tree_type="exon"
itree_edgelen$tree_type="intron"
stree_edgelen$tree_type="supercontig"


combined_edgelen=rbind(gtree_edgelen,etree_edgelen,itree_edgelen,stree_edgelen)
write.csv(combined_edgelen,file.path(stat_dir,"edgelengths.csv"),row.names=FALSE)

### ggarrange grid for node depth/branch length

bs.nd=full_nodekeys_df %>%
  mutate(tree_type = factor(tree_type, 
                            levels=c("genbank","exons","introns","supercontigs"),
                            labels=c("GenBank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=bootstrap))+
  geom_boxplot(aes(fill=tree_type),size=.25) +
  labs(x=NULL,y="Bootstrap (%)",fill="Gene subset") +
  scale_x_discrete(labels= c("deep","shallow")) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values=c("gray70","gray50","gray30","black"))

gcf.nd=full_nodekeys_df %>% filter(tree_type!="genbank") %>%
  #mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=gCF))+
  geom_boxplot(aes(fill=tree_type),size=.25,outlier.size=0.5) +
  labs(x="Node type",y="gCF (%)",fill="Gene subset") +
  scale_x_discrete(labels= c("deep","shallow")) +
  theme_classic(base_size=12) +
  scale_fill_manual(values=c("gray50","gray30","black"))


scf.nd=full_nodekeys_df %>% filter(tree_type!="genbank") %>%
  #mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
  ggplot(aes(x=parent_of_term,y=sCF))+
  geom_boxplot(aes(fill=tree_type),size=.25,outlier.size=0.5) +
  labs(x="Node type",y="sCF (%)",fill="Gene subset") +
  scale_x_discrete(labels= c("deep","shallow")) +
  theme_classic(base_size=12) +
  scale_fill_manual(values=c("gray50","gray30","black"))

bl=combined_edgelen %>% 
  mutate(tree_type = factor(tree_type, levels=c("GenBank","exons","introns","supercontigs"))) %>%
  ggplot()+
  geom_boxplot(aes(x=type,y=length,fill=tree_type),size=.25,outlier.size=0.5) +
  labs(x="Branch type",y="Length (subst/site)",fill="Dataset") +
  scale_fill_manual(values=c("gray70","gray50","gray30","black"))+
  theme_classic(base_size=12)


pdf(file=file.path(plot_dir,"nodeage_branchlen_box_grid_test2.pdf"),width=7.25,height=5.5)
fig8=ggpubr::ggarrange(bs.nd,gcf.nd,scf.nd,bl,labels=c("(a)","(b)","(c)","(d)"),label.y=0.13,
                  font.label=list(size=12),common.legend=TRUE) 
fig8
dev.off()
ggsave(filename="nodeage_branchlen_box_grid_test2.pdf",path=plot_dir,plot=fig8,device="pdf",
       width=7.25,height=5.5,units="in")

#######################
### non-figure stats ##
#######################
#### herbarium vs silica
hybpiper_stats=read.table(file.path(hyb_dir,"stats/hybpiper_stats_20200311.txt"),header=TRUE)
genbank_samples=unlist(read.csv(file.path(compare_dir,"genbank_compare_samples.txt"),header=FALSE))
libs_info=read.csv(file.path(hyb_dir,"libs_info_20200218.csv"),stringsAsFactors=FALSE)

hybpiper_stats_gsubset=dplyr::filter(hybpiper_stats,Name %in% genbank_samples)

hybpiper_stats_gs_libs=left_join(hybpiper_stats_gsubset,libs_info,by=c("Name"="Sample"))

se=function(x) sd(x)/sqrt(length(x))
table(hybpiper_stats_gs_libs$Source)
hybpiper_stats_gs_libs %>% group_by(Source) %>% summarize(mean.seqs=mean(GenesWithSeqs),se.seqs=se(GenesWithSeqs),mean.75=mean(GenesAt75pct),se.75=se(GenesAt75pct))

herb.seqs=hybpiper_stats_gs_libs %>% filter(Source=="herbarium") %>% dplyr::select(GenesWithSeqs,GenesAt75pct)
silica.seqs=hybpiper_stats_gs_libs %>% filter(Source=="silica") %>% dplyr::select(GenesWithSeqs,GenesAt75pct)

var.test(herb.seqs$GenesWithSeqs,silica.seqs$GenesWithSeqs) ## variance not significantly different, but t.test
t.test(herb.seqs$GenesWithSeqs,silica.seqs$GenesWithSeqs,var.equal = TRUE)
# Two Sample t-test
# 
# data:  herb.seqs$GenesWithSeqs and silica.seqs$GenesWithSeqs
# t = 0.67693, df = 77, p-value = 0.5005
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -5.20661 10.56979
# sample estimates:
#   mean of x mean of y 
# 307.6667  304.9851 

ggplot(data=hybpiper_stats_gs_libs) +
  geom_boxplot(aes(x=Source,y=GenesWithSeqs))

ggplot(data=hybpiper_stats_gs_libs) +
  geom_boxplot(aes(x=Source,y=GenesAt50pct))



### medians

all_nodekeys=group_by(all_nodekeys,dataset,tree_type)
medians=summarise(all_nodekeys,median.bs=median(bootstrap,na.rm = TRUE), median.gcf=median(gCF,na.rm=TRUE),median.scf=median(sCF,na.rm = TRUE))
range(medians$median.gcf,na.rm=TRUE)

quartiles.75=summarise(all_nodekeys,q75.bs=quantile(bootstrap,.75,na.rm = TRUE), q75.gcf=quantile(gCF,.75,na.rm=TRUE),q75.scf=median(sCF,.75,na.rm = TRUE))


# get max diff between filtering scheme medians for each gene subset and metric
medians %>% filter(tree_type!="genbank",dataset!="SortaDate_TL") %>% group_by(tree_type) %>% 
  summarise(diff.bs=max(median.bs)-min(median.bs),diff.gcf=max(median.gcf)-min(median.gcf),
            diff.scf=max(median.scf)-min(median.scf))

diff_subset=function(table,dat.col,comp.type){
  table=as.data.frame(table)
  #print(table[table$tree_type=="supercontigs",dat.col])
  #print(table[table$tree_type==comp.type,dat.col])
  return(table[table$tree_type=="supercontigs",dat.col]-table[table$tree_type==comp.type,dat.col])
}
# get max diff between gene subset medians for each gene filtering scheme and metric

diffs.bs=data.frame(dataset=character(0),sups_exons=numeric(0),sups_introns=numeric(0))
for(filt in c("full","intersection","SortaDate_BP","SortaDate_TL")){
  diff.exons=medians %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"median.bs","exons")
  print(diff.exons)
  diff.introns=medians %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"median.bs","introns")
  print(diff.introns)
  diffrow=data.frame(dataset=filt,sups_exons=diff.exons,sups_introns=diff.introns)
  diffs.bs=rbind(diffs.bs,diffrow)
}

diffs.gcf=data.frame(dataset=character(0),sups_exons=numeric(0),sups_introns=numeric(0))
for(filt in c("full","intersection","SortaDate_BP","SortaDate_TL")){
  print(filt)
  diff.exons=quartiles.75 %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"q75.gcf","exons")
  print(diff.exons)
  diff.introns=quartiles.75 %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"q75.gcf","introns")
  print(diff.introns)
  diffrow=data.frame(dataset=filt,sups_exons=diff.exons,sups_introns=diff.introns)
  diffs.gcf=rbind(diffs.gcf,diffrow)
}

diffs.scf=data.frame(dataset=character(0),sups_exons=numeric(0),sups_introns=numeric(0))
for(filt in c("full","intersection","SortaDate_BP","SortaDate_TL")){
  diff.exons=medians %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"median.scf","exons")
  print(diff.exons)
  diff.introns=medians %>% filter(tree_type!="genbank",dataset==filt) %>% group_by(dataset) %>% 
    diff_subset(.,"median.scf","introns")
  print(diff.introns)
  diffrow=data.frame(dataset=filt,sups_exons=diff.exons,sups_introns=diff.introns)
  diffs.scf=rbind(diffs.scf,diffrow)
}


nodetypes=all_nodekeys %>% group_by(parent_of_term,add = TRUE) %>% summarise(count.nodetype=n())

## look at enrichment stats (genbank subset)
## total bp of captured seqs
species_bp=read.csv(file.path(stat_dir,"species_bp_output.csv"),header=TRUE)
mean(species_bp$Exon.bp) # 147258.5
mean(species_bp$Intron.bp) # 221933

## genes (filter to genbank subset first)
gsubset_samples=read.table(file.path(compare_dir,"genbank_compare_samples.txt"))

summary_all=read.table(file.path(stat_dir,"veronica_summary_stats.txt"),header=TRUE)

summary_gsubset=dplyr::filter(summary_all,Name %in% gsubset_samples$V1)
dplyr::summarise_all(summary_gsubset[-1],mean)
# NumReads ReadsMapped PctOnTarget GenesMapped GenesWithContigs GenesWithSeqs
# 1  2797946    582314.4   0.2012785    335.6962              314           305
# GenesAt25pct GenesAt50pct GenesAt75pct Genesat150pct ParalogWarnings
# 1     272.7848     201.5823     113.1646     0.1265823        14.92405

#lengths_all=read.table(file.path(stat_dir,"veronica_seq_lengths.txt"),header=TRUE)
lengths_all=read.csv(file.path(stat_dir,"lengths_stats_20200311.csv"),header=TRUE)
rownames(lengths_all)=lengths_all[,1]
lengths_all_t=t(as.matrix(lengths_all[,-1]))
colnames(lengths_all_t)
zeros_50 = rownames(lengths_all_t)[which(lengths_all_t[,"percent_0s"]>.5)]
zeros_50=gsub("X","",zeros_50)
353-length(zeros_50)
length(which(lengths_all_t[,"percent_0s"]>=(1-(1/119))))

#The 79-species subset for the genbank comparison has 312(313?) genes that pass the 50% occupancy filter, and only 1 that doesn't make it from the 270 gene firstpass filter I used for the full dataset, and it's only one species short.


##### discordance factors (supplementary)

cf_tables=list()
for(name in full_tree_names){
  cf_fn=file.path(stat_dir,paste0("concordance_factors/",name,".cf.stat"))
  if(file.exists(cf_fn)){
    tree_type=gsub("Veronica_[0-9]+","",name)
    tree_type=gsub("_gsubset_concat","",tree_type)
    tree_type=gsub("_concat","",tree_type)
    print(tree_type)
    cf_table=read.table(cf_fn,header=TRUE)
    cf_table$tree_type=tree_type
    cf_tables[[tree_type]]=cf_table
  }
}

full_cf_tables=do.call("rbind",cf_tables)



png(filename = file.path(plot_dir,"gcf_gdf_full_box.png"),width=600,height=400)
full_cf_tables  %>% pivot_longer(cols = c(gCF,gDF1,gDF2,gDFP),
                                 names_to="gcf_type")  %>%
  group_by(tree_type,gcf_type) %>%
  ggplot(aes(x=gcf_type,y=value)) +
  geom_boxplot(aes(fill=tree_type)) +
  labs(x="Concordance factor",y="Percent of decisive loci",fill="Gene subset") +
  theme_classic(base_size=20) +
  scale_fill_manual(values=c("gray50","gray30","black"))
dev.off()

#####
# gdf test
cf_sup_fn=file.path(stat_dir,"concordance_factors/Veronica_229supercontigs_gsubset_concat.cf.stat")
cf_supercontigs=read.table(cf_sup_fn,header=TRUE)

boxplot(cf_supercontigs$gCF,cf_supercontigs$gDF1,cf_supercontigs$gDF2,cf_supercontigs$gDFP)
boxplot(cf_supercontigs$sCF,cf_supercontigs$sDF1,cf_supercontigs$sDF2)
boxplot(cf_supercontigs$gCF_N,cf_supercontigs$gN)
median(cf_supercontigs$gDFP)

cf_supercontigs[which(cf_supercontigs$gDF1>cf_supercontigs$gCF),]
cf_supercontigs[which(cf_supercontigs$gDF2>cf_supercontigs$gCF),]
#####

#################################
## low-paralog comparison plot ##
#################################
full_tree_names=c(supercontigs="Veronica_229supercontigs_gsubset_concat",
                  supercontigs_linsi_old="Veronica_gsubset_195supercontigs_noparalogs")
compare_nodekeys=compile_nodekeys(full_tree_names, "full",root=TRUE,plots=TRUE)
paralog_nodekeys_df=do.call("rbind",compare_nodekeys)

plot_dir=file.path(stat_dir,"manuscript_plots")

png(filename = file.path(plot_dir,"filter_compare_bootstrap_lowparalogs.png"),width=800,height=500)
bs=paralog_nodekeys_df %>%   
  mutate(tree_type = factor(tree_type, 
                            levels=c("supercontigs","supercontigs_noparalogs"),
                            labels=c("full","no_paralogs"))) %>%  
  ggplot(aes(x=tree_type,y=bootstrap)) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  labs(x=NULL,y="Bootstrap (%)") +
  theme_classic(base_size=25) +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())
#scale_fill_brewer(palette="Blues",direction=-1)
dev.off()

png(filename = file.path(plot_dir,"filter_compare_scf_lowparalogs.png"),width=600,height=400)
scf=paralog_nodekeys_df %>% 
  mutate(tree_type = factor(tree_type, 
                            levels=c("supercontigs","supercontigs_noparalogs"),
                            labels=c("full","no_paralogs"))) %>%  
  ggplot(aes(x=tree_type,y=sCF)) +
  geom_boxplot()  +
  labs(x="Supercontig data set",y="sCF (%)") +
  theme_classic(base_size=25) 
# scale_fill_brewer(palette="Blues",direction=-1)
dev.off()

png(filename = file.path(plot_dir,"filter_compare_gcf_lowparalogs.png"),width=800,height=500)
gcf=paralog_nodekeys_df   %>%
  mutate(tree_type = factor(tree_type, 
                            levels=c("supercontigs","supercontigs_noparalogs"),
                            labels=c("full","no_paralogs"))) %>%  
  ggplot(aes(x=tree_type,y=gCF)) +
  geom_boxplot() +
  labs(x="Supercontig data set",y="gCF (%)") +
  theme_classic(base_size=25) 
# theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
#scale_fill_brewer(palette="Blues",direction=-1)
dev.off()


png(filename = file.path(plot_dir,"compare_box_grid_lowparalogs_copyedit.png"),width=1000,height=750)
ggpubr::ggarrange(bs,gcf,scf,labels=c("(a)","(b)","(c)"),label.y=0.13,
                  font.label=list(size=20),legend="none") 
dev.off()



