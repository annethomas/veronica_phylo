library(dplyr)

getDescWrapper=function(tree,node,root){
  tree$tip.label[phytools::getDescendants(tree,node)[which(phytools::getDescendants(tree,node)<root)]]
}

## this is used in function, could add as argument
genbank_species=unlist(read.table(file.path(compare_dir,"genbank_compare_species.txt"),stringsAsFactors = FALSE))

add_node_group=function(tree,stat_table,info_table,group_name,group=NULL,subgroup=FALSE){
  ## checks and setup
  if(!is.null(group_name) & is.null(group)){
    #retrieve group species from info_table
    if(group_name %in% info_table$hebe_group){
      group=dplyr::filter(info_table, hebe_group==group_name & Species %in% genbank_species) %>% select("Species") %>% unlist()
    } else if(group_name %in% info_table$hebe_group_detail){
      group=dplyr::filter(info_table, hebe_group_detail==group_name & Species %in% genbank_species) %>% select("Species") %>% unlist()
    } else{
      print(unique(info_table$hebe_group))
      print(unique(info_table$hebe_group_detail))
      stop("please provide group_name present in info_table$hebe_group or info_table$hebe_group_detail")
      
    }
  } else if(is.null(group)){
    stop("Please provide group (vector) or group name")
  } else{
    if(any(!group %in% genbank_species)){
      print(group[which(!group %in% genbank_species)])
      stop("group includes species not in genbank_species")
    }
    print("using provided group species")
  }
  
  ## proceed
  print("species in target group:")
  print(group)
  mrca=getMRCA(tree,group)
  print(paste("mrca for",group_name,"is",mrca))
  
  if(length(getDescWrapper(tree,mrca,80))==length(group)){
    print("monophyletic")
    nodes=c(mrca,getDescendants(tree,mrca))
    print(nodes)
    
    # if(!all(c("node","group") %in% names(stat_table))){
    #   stop("please provide stat table with 'node' and 'group' columns")
    # }
    
    if(!"group" %in% names(stat_table)){
      #stop("please provide stat table with 'node' and 'group' columns")
      print("adding group column")
      stat_table$group=rep(NA,nrow(stat_table))
    }
    if(!"mrca" %in% names(stat_table)){
      #stop("please provide stat table with 'node' and 'group' columns")
      print("adding mrca column")
      stat_table$mrca=rep(FALSE,nrow(stat_table))
    }
    
    if(subgroup){
      if(!"subgroup" %in% names(stat_table)){
        #stop("please provide stat table with 'node' and 'group' columns")
        print("adding subgroup column")
        stat_table$subgroup=rep(NA,nrow(stat_table))
      }
      if(!"subgroup_mrca" %in% names(stat_table)){
        #stop("please provide stat table with 'node' and 'group' columns")
        print("adding subgroup_mrca column")
        stat_table$subgroup_mrca=rep(FALSE,nrow(stat_table))
      }
      stat_table[which(stat_table$node %in% nodes),"subgroup"]=group_name
      stat_table[which(stat_table$node == mrca),"subgroup_mrca"]=TRUE
    } else{
      stat_table[which(stat_table$node %in% nodes),"group"]=group_name
      stat_table[which(stat_table$node == mrca),"mrca"]=TRUE
    }
    
    return(stat_table)
    
  } else{ 
    print("paraphyletic species not in group:")
    print(setdiff(getDescWrapper(tree,mrca,80),group))
    stop("paraphyletic, please examine and break up group")
  }
}