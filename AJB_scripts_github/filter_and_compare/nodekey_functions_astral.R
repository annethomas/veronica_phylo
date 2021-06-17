library(ggplot2)
library(viridis)
library(dplyr)
library(ape)
library(phytools)


## parent of terminal node or deeper
## function to extract node information from rooted cf trees, print plots
make_nodekey_astral=function(cf_tree_rr,tree_name,plots=TRUE){
  genbank=FALSE
  astral=TRUE
  # same for genbank vs all
  nodekey=data.frame(node=1:cf_tree_rr$Nnode+Ntip(cf_tree_rr),nodelabel=cf_tree_rr$node.label,stringsAsFactors = FALSE)
  #print(paste("range of repeats of support combinations:",unique(table(nodekey$nodelabel)))) #1: no combinations are repeated
  
  cf_parentnodes=sapply(1:Ntip(cf_tree_rr),function(x) phytools::getParent(cf_tree_rr,x))
  #length(gcftree$tip.label) #alternative to ape::Ntip
  nodekey$parent_of_term=nodekey$node %in% cf_parentnodes
  # how many shallow and deep nodes
  print(table(nodekey$parent_of_term))

    nodekey=tidyr::separate(nodekey,col=nodelabel,into=c("posterior","gCF"),sep="/",remove=FALSE)
    nodekey$posterior=as.numeric(nodekey$posterior)
    nodekey$gCF=as.numeric(nodekey$gCF)


  
  if(plots){
    print("printing plots")
    # plot support values by node type/depth
    print(ggplot(nodekey, aes(x = parent_of_term, y = posterior)) + 
            ggtitle(tree_name) +
            geom_point() )
    print(ggplot(nodekey, aes(x = parent_of_term, y = posterior)) + 
            ggtitle(tree_name) +
            geom_boxplot() )
    
    
   
  }
  
  return(nodekey)
}

##################################
## root and plot with cf trees ##
##################################
compare_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/genbank_compare/"
stat_dir=file.path(compare_dir,"stats")
tree_astral_dir=file.path(compare_dir,"trees/astral/")
print(paste("warning: set compare_dir to",compare_dir,"and set stat_dir and tree_astral_dir"))

compile_nodekeys_astral=function(tree_names, dataset_type,plots=TRUE,root=TRUE,return_list=TRUE){
  ## read in trees with posterior and gcf values
  nodekey_list=list()
  tree_types=c("exons","introns","supercontigs")
  

  for(tree_name in tree_names){
    if(root){
      tree_fn=file.path(stat_dir,"concordance_factors/astral",paste0(tree_name,".cf.tree"))
      print(tree_fn)
      
      ## check file name (allow for genbank or other non-standard file)
      if(!file.exists(tree_fn)){
        
        print("Tree file doesn't exist, would you like to enter a different full filepath? Otherwise type no and this tree will be skipped.")
        tree_fn=readline("Reply: ")
        print(tree_fn)
        if(tree_fn=="no"){
          next
        } else if (!file.exists(tree_fn)){
          print("Tree file doesn't exist")
          next
        }
      }
      
      cf_tree=read.tree(tree_fn)
      
      cf_tree_rr = ape::root(cf_tree,edgelabel=TRUE,outgroup="Veronica_chamaedrys")
      plot(cf_tree_rr, type = "phylogram",cex=.7,show.node.label = TRUE)
      write.tree(cf_tree_rr,file=file.path(tree_astral_dir,paste0(tree_name,"_aperr.tre")))
      
    } else{
      tree_fn=file.path(tree_astral_dir,paste0(tree_name,"_aperr.tre"))
      print(tree_fn)
      
      ## check file name (allow for genbank or other non-standard file)
      if(!file.exists(tree_fn)){
        
        print("Tree file doesn't exist, would you like to enter a different full filepath? Otherwise type no and this tree will be skipped.")
        tree_fn=readline("Reply: ")
        print(tree_fn)
        if(tree_fn=="no"){
          next
        } else if (!file.exists(tree_fn)){
          print("Tree file doesn't exist")
          next
        }
      }
      
      cf_tree_rr=read.tree(tree_fn)
      plot(cf_tree_rr, type = "phylogram",cex=.7,show.node.label = TRUE)
      
    }
    
    nodekey=make_nodekey_astral(cf_tree_rr,tree_name,plots=plots)
    nodekey$dataset=dataset_type
    tree_type=names(which(sapply(tree_types,grepl, tree_name)))
    if(length(tree_type)!=1){
      print(tree_type)
      print("Don't know tree type, would you like to enter it? Otherwise type no and this tree will be skipped.")
      tree_type=readline("Reply: ")
      if(tree_type=="no"){
        next
      } 
        #print("Don't know tree type, skipping") ##could add interactive bit...
        #next
    }
    nodekey$tree_type=tree_type
    nodekey_list[[tree_name]]=nodekey
    

  }
  
  ## rbind
  all_nodekeys=do.call("rbind",nodekey_list)

  print(
    all_nodekeys  %>%
      # mutate(tree_type = factor(tree_type, levels=c("genbank","exons","introns","supercontigs"))) %>%
      ggplot(aes(x=tree_type,y=posterior)) +
      geom_boxplot()
    #  geom_boxplot(aes(fill=tree)) #+
    #  theme(axis.text.x = element_text(angle = 90))
    ## specify colors to match the three-subset figs
    ## add count as text annotation
  )
  
  print(  
    all_nodekeys %>% 
      ggplot(aes(x=tree_type,y=gCF))+
      geom_boxplot()
  )

  
  if(return_list){
    return(nodekey_list)
  }else{
    return(all_nodekeys)
  }
}
