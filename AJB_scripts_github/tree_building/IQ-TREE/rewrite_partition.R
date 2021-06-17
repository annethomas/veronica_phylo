compare_dir="~/HybPiper/genbank_compare/"

tree_dir=file.path(compare_dir,"iqtree/")

#gene_pos_fn=commandArgs(trailingOnly=TRUE)[1] ## eg "Veronica_257exons_concat_gene_pos.csv"
#old_partition_fn=commandArgs(trailingOnly=TRUE)[2]
#new_partition_fn=commandArgs(trailingOnly=TRUE)[3]

old_basename=commandArgs(trailingOnly=TRUE)[1] ## eg Veronica_177exons_gsubset_concat
new_basename=commandArgs(trailingOnly=TRUE)[2] ## eg Veronica_175exons_gsubset_concat

old_concat_dir=file.path(tree_dir,old_basename)
new_concat_dir=file.path(tree_dir,new_basename)

old_partition=file(file.path(old_concat_dir,paste0(old_basename,"_gene_pos.nex.best_scheme.nex")),"r")
new_partition=file(file.path(new_concat_dir,paste0(new_basename,".best_scheme.nex")),"a")

gene_pos_fn=file.path(new_concat_dir,paste0(new_basename,"_gene_pos.csv"))
gene_pos = read.csv(gene_pos_fn,stringsAsFactors = FALSE)
genelist=gene_pos$gene



while ( TRUE ) {
  line = readLines(old_partition, n = 1)
  if ( length(line) == 0 ) {
    break
  } else if(substr(line,3,nchar("  charset"))=="charset"){
    
    print(line)
    genes=stringr::str_extract(substr(line,nchar("  charset "),nchar(line)),"[gene,0-9,_]+")
    gene1 = stringr::str_extract(genes,"[0-9]+")
    if(as.numeric(gene1) %in% genelist){
      #keep this partition and extract and rewrite gene positions
      gene_set = strsplit(genes,split="_")
      pos_list=c()
      for(g in gene_set[[1]]){
        print(g)
        gene = as.numeric(substr(g,5,8))
        pos_list=c(pos_list,paste0(gene_pos[gene_pos$gene==gene,"start"],"-",gene_pos[gene_pos$gene==gene,"end"]))
      }
      print(pos_list)
      #append first part of line and pos_list to new file, plus ;
      keep = stringr::str_extract(line,".+=")
      new_line = paste0(keep," ",paste(pos_list,collapse = "  "),";")
      writeLines(new_line,new_partition)
    }  else{print("not included")}
    
  } else{
    ## for checking and removing rejected genes from model assignments
    print(line)
    genes = stringr::str_extract(line,": [gene,0-9,_]+")
    
    if(!is.na(genes)){
      gene1 = stringr::str_extract(genes,"[0-9]+")
      if(as.numeric(gene1) %in% genelist){
        writeLines(line,new_partition)
      } else{
        print("not included")
        ## check that last partition isn't rejected; otherwise need to add in final semicolon
      }
    } else{
      #it's a formatting line, should be included
      writeLines(line,new_partition)
    }
  }
  
}

close(old_partition)
close(new_partition)


