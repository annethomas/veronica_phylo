
exon_tree=ape::read.tree(file.path(iq_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat_allbstrees.suptree")) 
genbank_tree=ape::read.tree(file.path(iq_dir,"Veronica_genbank_concat5/Veronica_genbank_concat5_allbstrees.suptree"))
library(TreeDist)
#https://ms609.github.io/TreeDist/reference/index.html#section-tree-distance-measures

hyb_dir="~/HybPiper/genbank_compare/"
iq_dir=file.path(hyb_dir,"iqtree")

exon_tree$tip.label[which(exon_tree$tip.label=="Veronica_hookeriana-BDM12")]="Veronica_hookeriana"
exon_tree$tip.label[which(exon_tree$tip.label=="Veronica_vernicosa-2926")]="Veronica_vernicosa"

####################
##Distance metrics##
####################
## doesn't matter what order the trees are in

### Information/entropy-based
TreeDistance(genbank_tree,exon_tree) # 0.4598027 
MutualClusteringInfo(exon_tree,genbank_tree,normalize=TRUE) # 0.4598027
ClusteringInfoDistance(exon_tree,genbank_tree,normalize=TRUE) # 0.5401973

PhylogeneticInfoDistance(exon_tree,genbank_tree,normalize=TRUE) # 0.6021813
DifferentPhylogeneticInfo(exon_tree,genbank_tree, normalize=TRUE) # 0.6021813
SharedPhylogeneticInfo(exon_tree,genbank_tree,normalize=TRUE) # 0.3978187

MatchingSplitInfo(exon_tree,genbank_tree,normalize=TRUE) # 0.6856767
MatchingSplitInfoDistance(exon_tree,genbank_tree,normalize=TRUE) # 0.3143233

## non-normalized versions
MutualClusteringInfo(exon_tree,genbank_tree,normalize=FALSE) # 13.95318
ClusteringInfoDistance(exon_tree,genbank_tree,normalize=FALSE) # 32.78568

PhylogeneticInfoDistance(exon_tree,genbank_tree,normalize=FALSE) # 2445.644
DifferentPhylogeneticInfo(exon_tree,genbank_tree, normalize=FALSE) # 2445.644
SharedPhylogeneticInfo(exon_tree,genbank_tree,normalize=FALSE) # 807.8324

MatchingSplitInfo(exon_tree,genbank_tree,normalize=FALSE) # 1392.373
MatchingSplitInfoDistance(exon_tree,genbank_tree,normalize=FALSE) # 1276.564

## comparison to random trees:
ExpectedVariation(exon_tree, genbank_tree, samples = 10000)
#                             Estimate   Std. Err.         sd     n
# SharedPhylogeneticInfo      93.53897 0.123638562 12.3638562 10000
# MatchingSplitInfo          854.08036 0.211381149 21.1381149 10000
# MutualClusteringInfo         3.18964 0.002491386  0.2491386 10000
# DifferentPhylogeneticInfo 3874.23122 0.247277124 24.7277124 10000
# MatchingSplitInfoDistance 2353.14843 0.422762298 42.2762298 10000
# ClusteringInfoDistance      54.31276 0.004982771  0.4982771 10000

### Generalized Robinson-Foulds
NyeSimilarity(exon_tree, genbank_tree,similarity=FALSE, normalize=TRUE) # 0.4856765
NyeSimilarity(exon_tree, genbank_tree,similarity=TRUE, normalize=TRUE) # 0.5143235

NyeSimilarity(exon_tree, genbank_tree,similarity=FALSE, normalize=FALSE) # 73.82283
NyeSimilarity(exon_tree, genbank_tree,similarity=TRUE, normalize=FALSE) # 39.08859

JaccardRobinsonFoulds(exon_tree,genbank_tree,normalize=TRUE) # 0.4856765
JaccardRobinsonFoulds(exon_tree,genbank_tree,normalize=FALSE) # 73.82283

RobinsonFoulds(exon_tree,genbank_tree,normalize=TRUE) # 0.8026316
RobinsonFoulds(exon_tree,genbank_tree,normalize=FALSE) # 122
InfoRobinsonFoulds(exon_tree,genbank_tree,normalize=TRUE) # 0.8485187
InfoRobinsonFoulds(exon_tree,genbank_tree,normalize=FALSE) # 3446.097


MatchingSplitDistance(exon_tree,genbank_tree,normalize=FALSE) # 415

KendallColijn(exon_tree,genbank_tree) # 442.5912

MASTSize(exon_tree, genbank_tree,rooted=FALSE) # 30
MASTInfo(exon_tree, genbank_tree,rooted=FALSE) # 122.7083

NNIDist(exon_tree, genbank_tree)
#      lower  best_lower tight_upper  best_upper loose_upper  fack_upper	li_upper
#         64          81          NA         325         356         368     369

SPRDist(exon_tree,genbank_tree) # 27

PathDist(exon_tree, genbank_tree) # 357.7793
## "Use of the path distance is discouraged as it emphasizes shallow relationships at the expense of deeper (and arguably more fundamental) relationships (Farris, 1973)."

######################
##adding more trees ##
######################
exon_tree=ape::read.tree(file.path(iq_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat_allbstrees.suptree")) 
genbank_tree=ape::read.tree(file.path(iq_dir,"Veronica_genbank_concat5/Veronica_genbank_concat5_allbstrees.suptree"))

intron_tree=ape::read.tree(file.path(iq_dir,"Veronica_167introns_gsubset_concat/Veronica_167introns_gsubset_concat_allbstrees.suptree"))

intron_tree$tip.label[which(intron_tree$tip.label=="Veronica_hookeriana-BDM12")]="Veronica_hookeriana"
intron_tree$tip.label[which(intron_tree$tip.label=="Veronica_vernicosa-2926")]="Veronica_vernicosa"

supercontig_tree=ape::read.tree(file.path(iq_dir,"Veronica_229supercontigs_gsubset_concat/Veronica_229supercontigs_gsubset_concat_allbstrees.suptree"))



TreeDistance(list(g=genbank_tree,e=exon_tree,i=intron_tree,s=supercontig_tree))
#           g         e         i         s
# g 1.0000000 0.4598027 0.4683749 0.4708813
# e 0.4598027 1.0000000 0.5539304 0.5989622
# i 0.4683749 0.5539304 1.0000000 0.6327310
# s 0.4708813 0.5989622 0.6327310 1.0000000

MutualClusteringInfo(list(g=genbank_tree,e=exon_tree,i=intron_tree,s=supercontig_tree))
#          [,1]     [,2]     [,3]     [,4]
# [1,] 33.44465 13.95318 14.30576 14.13307
# [2,] 13.95318 27.24739 15.20249 16.12134
# [3,] 14.30576 15.20249 27.64214 17.15512
# [4,] 14.13307 16.12134 17.15512 26.58351

ExpectedVariation(exon_tree,intron_tree,samples=10000)

#                             Estimate   Std. Err.         sd     n
#SharedPhylogeneticInfo      86.691017 0.122435788 12.2435788 10000
#MatchingSplitInfo          747.465997 0.177503774 17.7503774 10000
#MutualClusteringInfo         3.355208 0.002398102  0.2398102 10000
#DifferentPhylogeneticInfo 3400.807902 0.244871577 24.4871577 10000
#MatchingSplitInfoDistance 2079.257943 0.355007547 35.5007547 10000
#ClusteringInfoDistance      48.179114 0.004796203  0.4796203 10000

##comparing different versions of same dataset
old_intron_tree=ape::read.tree(file.path(iq_dir,"Veronica_169introns_gsubset_concat/Veronica_169introns_gsubset_concat_allbstrees.suptree"))

old_intron_tree$tip.label[which(old_intron_tree$tip.label=="Veronica_hookeriana-BDM12")]="Veronica_hookeriana"
old_intron_tree$tip.label[which(old_intron_tree$tip.label=="Veronica_vernicosa-2926")]="Veronica_vernicosa"

TreeDistance(intron_tree,old_intron_tree)
# 0.7017369
MutualClusteringInfo(intron_tree,old_intron_tree)
# 19.20321



## comparing intersection dataset
library(TreeDist)
hyb_dir="~/HybPiper/genbank_compare/"
iq_dir=file.path(hyb_dir,"iqtree")

genbank_tree=ape::read.tree(file.path(iq_dir,"Veronica_genbank_concat5/Veronica_genbank_concat5_allbstrees.suptree"))
exon_tree=ape::read.tree(file.path(iq_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat_allbstrees.suptree")) 
intron_tree=ape::read.tree(file.path(iq_dir,"Veronica_167introns_gsubset_concat/Veronica_167introns_gsubset_concat_allbstrees.suptree"))
supercontig_tree=ape::read.tree(file.path(iq_dir,"Veronica_229supercontigs_gsubset_concat/Veronica_229supercontigs_gsubset_concat_allbstrees.suptree"))

## intersection dataset
etree=ape::read.tree(file.path(iq_dir,"Veronica_exons_intersect_concat/Veronica_exons_intersect_concat_allbstrees.suptree"))
itree=ape::read.tree(file.path(iq_dir,"Veronica_introns_intersect_concat/Veronica_introns_intersect_concat_allbstrees.suptree"))
stree=ape::read.tree(file.path(iq_dir,"Veronica_supercontigs_intersect_concat/Veronica_supercontigs_intersect_concat_allbstrees.suptree"))

## sortadate "strict" dataset (median cutoff for bipartition support and tree length) [number of genes in file names]
esstree=ape::read.tree(file.path(iq_dir,"Veronica_med51exons_gsubset_concat/Veronica_med51exons_gsubset_concat_allbstrees.suptree"))
isstree=ape::read.tree(file.path(iq_dir,"Veronica_med63introns_gsubset_concat/Veronica_med63introns_gsubset_concat_allbstrees.suptree"))
ssstree=ape::read.tree(file.path(iq_dir,"Veronica_med80supercontigs_gsubset_concat/Veronica_med80supercontigs_gsubset_concat_allbstrees.suptree"))

## sortadate "lax" datset (.02 bipartition dataset (roughly median for all), .4 tree length (none to low number removed))
esltree=ape::read.tree(file.path(iq_dir,"Veronica_fix96exons_gsubset_concat/Veronica_fix96exons_gsubset_concat_allbstrees.suptree"))
isltree=ape::read.tree(file.path(iq_dir,"Veronica_fix116introns_gsubset_concat/Veronica_fix116introns_gsubset_concat_allbstrees.suptree"))
ssltree=ape::read.tree(file.path(iq_dir,"Veronica_fix150supercontigs_gsubset_concat/Veronica_fix150supercontigs_gsubset_concat_allbstrees.suptree"))


TreeDistance(list(g=genbank_tree,e.full=exon_tree,i.full=intron_tree,s.full=supercontig_tree,e.int=etree,i.int=itree,s.int=stree))
TreeDistance(list(g=genbank_tree,e.full=exon_tree,i.full=intron_tree,s.full=supercontig_tree,e.int=etree,i.int=itree,s.int=stree,
			e.sorts=esstree,i.sorts=isstree,s.sorts=ssstree,e.sortl=esltree,i.sortl=isltree,s.sortl=ssltree))



rooted=ape::read.tree("/home/at820/HybPiper/genbank_compare/stats/concordance_factors/Veronica_229supercontigs_aperr.tre")
unrooted=ape::read.tree("/home/at820/HybPiper/genbank_compare/stats/concordance_factors/Veronica_229supercontigs_gsubset_concat_test.cf.tree")
TreeDistance(unrooted, rooted)
#[1] 1



### 202011 New supercontig all-taxa tree
library(TreeDist)
#https://ms609.github.io/TreeDist/reference/index.html#section-tree-distance-measures

hyb_dir="~/HybPiper/"
iq_dir=file.path(hyb_dir,"iqtree")

tree1=ape::read.tree(file.path(iq_dir,"Veronica_239supercontigs/Veronica_239supercontigs_gene_pos.nex.treefile")) 
tree2=ape::read.tree(file.path(iq_dir,"Veronica_239supercontigs/Veronica_239supercontigs_MLforBS.treefile"))


TreeDistance(tree1,tree2) #0.64

## 20210113 ASTRAL version of filtering schemes
library(TreeDist)
#https://ms609.github.io/TreeDist/reference/index.html#section-tree-distance-measures
hyb_dir="~/HybPiper/"


compare_dir="~/HybPiper/genbank_compare/"
astral_dir=file.path(compare_dir,"astral")
iq_dir=file.path(compare_dir,"iqtree")

genbank_tree=ape::read.tree(file.path(iq_dir,"Veronica_genbank_concat5/Veronica_genbank_concat5_allbstrees.suptree"))
exon_tree=ape::read.tree(file.path(astral_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat.astral.tre")) 
intron_tree=ape::read.tree(file.path(astral_dir,"Veronica_167introns_gsubset_concat/Veronica_167introns_gsubset_concat.astral.tre"))
supercontig_tree=ape::read.tree(file.path(astral_dir,"Veronica_229supercontigs_gsubset/Veronica_229supercontigs_gsubset.astral.tre"))

## intersection dataset
etree=ape::read.tree(file.path(astral_dir,"Veronica_exons_intersect_concat/Veronica_exons_intersect_concat.astral.tre"))
itree=ape::read.tree(file.path(astral_dir,"Veronica_introns_intersect_concat/Veronica_introns_intersect_concat.astral.tre"))
stree=ape::read.tree(file.path(astral_dir,"Veronica_supercontigs_intersect_concat/Veronica_supercontigs_intersect_concat.astral.tre"))

## sortadate "strict" dataset (median cutoff for bipartition support and tree length) [number of genes in file names]
esstree=ape::read.tree(file.path(astral_dir,"Veronica_med51exons_gsubset_concat/Veronica_med51exons_gsubset_concat.astral.tre"))
isstree=ape::read.tree(file.path(astral_dir,"Veronica_med63introns_gsubset_concat/Veronica_med63introns_gsubset_concat.astral.tre"))
ssstree=ape::read.tree(file.path(astral_dir,"Veronica_med80supercontigs_gsubset_concat/Veronica_med80supercontigs_gsubset_concat.astral.tre"))

## sortadate "lax" datset (.02 bipartition dataset (roughly median for all), .4 tree length (none to low number removed))
esltree=ape::read.tree(file.path(astral_dir,"Veronica_fix96exons_gsubset_concat/Veronica_fix96exons_gsubset_concat.astral.tre"))
isltree=ape::read.tree(file.path(astral_dir,"Veronica_fix116introns_gsubset_concat/Veronica_fix116introns_gsubset_concat.astral.tre"))
ssltree=ape::read.tree(file.path(astral_dir,"Veronica_fix150supercontigs_gsubset_concat/Veronica_fix150supercontigs_gsubset_concat.astral.tre"))


TreeDistance(list(g=genbank_tree,e.full=exon_tree,i.full=intron_tree,s.full=supercontig_tree,e.int=etree,i.int=itree,s.int=stree))
TreeDistance(list(g=genbank_tree,e.full=exon_tree,i.full=intron_tree,s.full=supercontig_tree,e.int=etree,i.int=itree,s.int=stree,
			e.sorts=esstree,i.sorts=isstree,s.sorts=ssstree,e.sortl=esltree,i.sortl=isltree,s.sortl=ssltree))


### Astral vs iqtree

g=ape::read.tree(file.path(iq_dir,"Veronica_genbank_concat5/Veronica_genbank_concat5_allbstrees.suptree"))
e.iq=ape::read.tree(file.path(iq_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat_allbstrees.suptree")) 
i.iq=ape::read.tree(file.path(iq_dir,"Veronica_167introns_gsubset_concat/Veronica_167introns_gsubset_concat_allbstrees.suptree"))
s.iq=ape::read.tree(file.path(iq_dir,"Veronica_229supercontigs_gsubset_concat/Veronica_229supercontigs_gsubset_concat_allbstrees.suptree"))

e.ast=ape::read.tree(file.path(astral_dir,"Veronica_173exons_gsubset_concat/Veronica_173exons_gsubset_concat.astral.tre")) 
i.ast=ape::read.tree(file.path(astral_dir,"Veronica_167introns_gsubset_concat/Veronica_167introns_gsubset_concat.astral.tre"))
s.ast=ape::read.tree(file.path(astral_dir,"Veronica_229supercontigs_gsubset/Veronica_229supercontigs_gsubset.astral.tre"))

TreeDistance(list(g=g,e.iq=e.iq,i.iq=i.iq,s.iq=s.iq,e.ast=e.ast,i.ast=i.ast,s.ast=s.ast))




