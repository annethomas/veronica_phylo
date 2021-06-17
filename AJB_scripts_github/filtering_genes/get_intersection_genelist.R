compare_dir="~/HybPiper/genbank_compare/"
data_dir="/production/at820/genbank_subset"

exons173=unlist(read.table(file.path(data_dir,"exons_gsubset_trm_cleaned_tiptrim_secondpass3.txt")))
supercontigs229=unlist(read.table(file.path(data_dir,"supercontigs_gsubset_trm_cleaned_tiptrim_secondpass.txt")))
introns167=unlist(read.table(file.path(data_dir,"introns_gsubset_trm_cleaned_tiptrim_secondpass.txt")))


print(length(intersect(intersect(exons173,introns167),supercontigs229))) #119 (51% of full supercontig list)

intersect_genes=intersect(intersect(exons173,introns167),supercontigs229)
write.table(intersect_genes,file.path(data_dir,"secondpass_intersection_genes.txt"),row.names=FALSE,col.names=FALSE)
