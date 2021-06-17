library(ape)

hyb_dir="~/HybPiper"
iq_dir=file.path(hyb_dir,"iqtree")

tree_base=commandArgs(trailingOnly=TRUE)[1] ## e.g. Veronica_239supercontigs
tree_mod=commandArgs(trailingOnly=TRUE)[2] ## everything after base but before .treefile, e.g. "_MLforBS" (will need to modify if using trees not ending in treefile)
tree_suffix=commandArgs(trailingOnly=TRUE)[3] ## e.g. .suptree or .treefile
tree_dir=commandArgs(trailingOnly=TRUE)[4] ## full path of input dir, if different from iqtree/tree_base; otherwise can leave out

if(is.na(tree_dir)){
	tree_dir=file.path(iq_dir,tree_base)
}


tree_file=file.path(tree_dir,paste0(tree_base,tree_mod,tree_suffix))
print(tree_file)

r_tree_file=file.path(tree_dir,paste0(tree_base,tree_mod,"_rooted",tree_suffix))
if(file.exists(r_tree_file)){
	stop(paste(r_tree_file,"already exists"))
} else{
	tree=read.tree(tree_file)
	tree.r=root(tree,"S4D_V_chamaedrys",edgelabel=TRUE,resolve.root=TRUE)

	print(r_tree_file)
	write.tree(tree.r,r_tree_file)
}

