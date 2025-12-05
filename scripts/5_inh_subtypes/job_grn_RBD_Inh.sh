#!/bin/bash
#SBATCH --array 0
#SBATCH --job-name grn_RBD
#SBATCH --mem 100g -t 12:00:00
#SBATCH --partition=scavenge
#SBATCH --cpus-per-task 5


# DO NOT EDIT LINE BELOW
module load miniconda
conda activate conda_envs/scglue
python ../../helpers/save_grn.py \
  --loom ../../Signac_plots/saved_files/rbd_inh_all_tfs_rna_INH_SUBSET.loom \
  --tfs ../../Signac_plots/saved_files/all_tfs.txt \
  --draft ../../data/GRN_Inh_subset/rbd_draft_grn.csv \
  --links ../../Signac_plots/saved_files/RBD_gene_peak_links_inh_subset.tsv \
  --peaks ../../Signac_plots/saved_files/RBD_peaks_inh_subset.bed \
  --rankings ../../data/GRN_Inh_subset/rbd_signac.genes_vs_tracks.rankings.feather \
  --rankings_prom ../../data/GRN_Inh_subset/rbd_supp.genes_vs_tracks.rankings.feather \
  --annotation ../../data/GRN_Inh_subset/rbd_ctx_annotation.tsv \
  --pruned ../../data/GRN_Inh_subset/RBD_pruned_grn.csv \
  --gene2peak-out ../../data/GRN_Inh_subset/RBD_gene2peak_edges.csv \
  --peak2tf-out ../../data/GRN_Inh_subset/RBD_peak2tf_edges.csv \
  --flank2tf-out ../../data/GRN_Inh_subset/RBD_flank2tf_edges.csv