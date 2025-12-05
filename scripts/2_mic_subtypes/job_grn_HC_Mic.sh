#!/bin/bash
#SBATCH --array 0
#SBATCH --job-name grn_HC
#SBATCH --mem 100g -t 6:00:00
#SBATCH --partition=scavenge
#SBATCH --cpus-per-task 5


# DO NOT EDIT LINE BELOW
module load miniconda
conda activate conda_envs/scglue
python ../../helpers/save_grn.py \
  --loom ../../Signac_plots/saved_files/hc_mic_all_tfs_rna_MIC.loom \
  --tfs ../../Signac_plots/saved_files/all_tfs.txt \
  --draft ../../data/GRN_Mic/hc_draft_grn.csv \
  --links ../../Signac_plots/saved_files/HC_gene_peak_links_mic.tsv \
  --peaks ../../Signac_plots/saved_files/HC_peaks_mic.bed \
  --rankings ../../data/GRN_Mic/hc_signac.genes_vs_tracks.rankings.feather \
  --rankings_prom ../../data/GRN_Mic/hc_supp.genes_vs_tracks.rankings.feather \
  --annotation ../../data/GRN_Mic/hc_ctx_annotation.tsv \
  --pruned ../../data/GRN_Mic/HC_pruned_grn.csv \
  --gene2peak-out ../../data/GRN_Mic/HC_gene2peak_edges.csv \
  --peak2tf-out ../../data/GRN_Mic/HC_peak2tf_edges.csv \
  --flank2tf-out ../../data/GRN_Mic/HC_flank2tf_edges.csv