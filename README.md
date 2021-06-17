# veronica_phylo
Scripts for work on the phylogeny of New Zealand Veronica based on Angiosperms353 target enrichment data

This is the rough outline of the order of scripts in the pipeline for the analyses in Thomas et al 2021 AJB. The scripts were tailored to a university cluster (e.g. job submission with a script not included) and the pipeline was not entirely automated due to the necessity of running non-synchronous jobs.

## Scraping Genbank for existing Veronica sequences and preparing alignments
### script directory: Genbank_sequences_pipeline
_Not included yet: run pyPhlawd_
#### align_clusters.sh
#### run_trimal_clusters.sh
#### concat_supermatrix_write_partition_clusters.R
#### run_modelfinder_cluster.sh
#### get_genbank_accessions.sh

## Trimming and quality checking raw Illumina sequences available in the Sequence Read Archive
### script directory: read_processing
#### trimmomatic_batch.sh
#### fastqc_batch.sh

## Wrappers for running HybPiper
### script directory: HybPiper_wrappers
#### hybpiper_batch_readsfirst_all.sh
submits job for each specimen to run hybpiper (cluster-specific submission script not included)
#### intronerate_batch.sh
loops through specimens with intronerate.py
#### paralog_investigate_loop.sh
saves output of paralog_investigator.py
#### run_retrieve_sequences.sh
retrieves sequences for exon nucleotides, exon amino acids, introns, and supercontigs
#### get_hypbiper_stats.sh
runs HybPiper functions for collating stats on retrieved genes

## Aligning, trimming, and cleaning Angiosperms353 sequences
### script dir: align_and_clean
#### create_firstpass_genelist.R
filters unaligned gene files (outputting list of genes) based on HybPiper stats of recovery and paralog warnings
#### run_mafft.sh <type>
aligns sequences for each gene given the type of sequences (exons_FNA, exons_FAA, introns, supercontig)
#### run_pal2nal_trimal.sh
completes in-frame alignment of exon nucleotides using exon amino acids; trims alignments
#### run_trimal.sh <type>
trims alignments for each gene given the type of sequences (introns, supercontig)
#### remove_short_sequences.R <input_dir> <filter_file>
removes stub sequences of less than 5% of alignment length given alignment directory and list of genes; writes to new directory

## Trimming long tips from gene trees and updating alignments (a step in alignment cleaning)
### script dir: align_and_clean
#### run_iqtree_genetrees_loop_wrapper.sh
builds a tree for each gene with IQTREE
#### _Manually create trimmed_trees dir inside of above genetrees/<filebase> and cp *.treefile from the iqtree output into trimmed_trees_
#### trim_tips_edit.py <genetree_dir> <file_ending> .2 .5
edited version of Yang and Smith pipeline script for removing long tips; this version also removes long internal branches and subtending tips
#### trim_tips_for_aln.R <trim_tree_dir> <tree_suffix> <aln_dir> <aln_suffix> <out_aln_dir> <outgroups_to_keep…>
removes sequences of trimmed tips from alignments and save to new directory (out_aln_dir)
#### _Rerun run_iqtree_genetrees_loop_wrapper.sh with new alignments_

## Alignment-stats-based filtering ("second pass")
### script dir: align_and_clean
#### create_secondpass_genelist_wrapper.R <input_dir>
generates alignment statistics and creates list of filtered alignments based on specified thresholds (creates the baseline dataset for later filtering); calls create_secondpass_genelist.R and get_gene_stats.R

## Tree building
### script dir: tree_building/IQ-TREE (note: replace path to IQTREE installation in scripts)
#### rename_tips.R
rename tip labels in alignment or tree from sample IDs to species names
#### concat_supermatrix_write_partition.R
combine individual gene alignments to supermatrix and keep track of positions to write gene partition file
#### run_modelfinder.sh
merge partitions and find best model for each (limited to 10% of search space)
#### run_iqtree_bootstrap_loop.sh
run 100 bootstraps
#### rewrite_partition.sh
remove short alignments remaining after merging partitions that sometimes prevent bootstraps from converging
#### run_MLforBS.sh
run an independent ML tree to map bootstrap values on
#### prep_bootstrap_trees.sh
concatenate bootstrap trees, calculate percentages and annotate ML tree
### script dir: tree_building/ASTRAL_SVDQUARTETS
#### run_astral_all.sh
run astral and gcf for all filtering schemes
#### run_astral.sh
single astral run with automatic prep steps
#### run_paup_svdquartets.sh
very basic svdquartets job submission
#### run_iqtree_gcf.sh
calculate gCF, sCF for a particular species tree given gene trees or alignments

## Filtering genes (additional filtering schemes after "second pass" above) and comparing trees; final plots
###	script dir: filter_and_compare
#### get_intersection_genelist.R
determine which genes are present in all gene sets (exons, introns, supercontigs) after second pass filtering for intersection dataset
#### run_sortadate.sh
Run SortaDate to get list of metrics for filtering gene trees
#### create_sortadate_filtering_schemes.R
use sortadate stats to make gene lists for filtering schemes (dataset names may not reflect manuscript versions); repeat tree building steps above with new gene lists
#### run_TreeDist.R
run various tree distance metrics including MutualClusteringInfo
#### 2021AJB_manuscript_plots.R
prepare final data/tree comparisons and create all R-generated figures for manuscript as well as additional stats (sources nodekey_functions.R and subtree_functions.R)
