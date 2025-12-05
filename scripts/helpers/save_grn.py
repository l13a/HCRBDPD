#!/usr/bin/env python3
"""
run_signac_scenic.py

Construct a SCENIC-style TF→gene regulatory network using Signac-derived gene-peak links
and pyscenic, with all file paths supplied via command-line arguments.
"""
import argparse
import subprocess
import pandas as pd
import networkx as nx
import numpy as np
import scglue
import loompy
import scanpy as sc 


def main():
    parser = argparse.ArgumentParser(
        description="Build draft GRN, rank TF-gene links via Signac peaks, and prune with pyscenic"
    )
    parser.add_argument("--loom",        help="Input loom file (expression)",               required=True)
    parser.add_argument("--tfs",         help="Input TF list txt file",                   required=True)
    parser.add_argument("--draft",       help="Output draft GRN CSV file",                required=True)
    parser.add_argument("--links",       help="Input gene-peak links TSV file",           required=True)
    parser.add_argument("--peaks",       help="Input peaks BED file",                     required=True)
    parser.add_argument("--rankings",    help="Output cis-regulatory rankings feather file", required=True)
    parser.add_argument("--rankings_prom",    help="Output promoter rankings feather file", required=True)
    parser.add_argument("--annotation",  help="Output ctx annotation TSV file",           required=True)
    parser.add_argument("--pruned",      help="Output pruned GRN CSV file",               required=True)
    parser.add_argument("--motif-bed",   help="Input JASPAR/ChIP BED file",              default="bed_files/JASPAR2022-hg38.bed.gz")
    parser.add_argument("--seed",        type=int, help="Random seed for reproducibility",        default=0)
    parser.add_argument("--workers",     type=int, help="Number of workers for pyscenic",       default=20)
    parser.add_argument("--gene2peak-out", help="Save gene–peak edge list CSV", required=False)
    parser.add_argument("--peak2tf-out",   help="Save peak–TF edge list CSV",   required=False)
    parser.add_argument("--flank2tf-out",   help="Save flank–TF edge list CSV",   required=False)
    args = parser.parse_args()

    # 1. Draft network: GRNBoost2 via pyscenic
    cmd_grn = [
        "pyscenic", "grn",
        args.loom, args.tfs,
        "-o", args.draft,
        "--seed", str(args.seed),
        "--num_workers", str(args.workers),
        "--cell_id_attribute", "CellID",
        "--gene_attribute", "Gene"
    ]
    print("[1/6] Running GRNBoost2 draft network...")
    subprocess.run(cmd_grn, check=True)

    # 2. Build gene–peak graph from Signac links
    print("[2/6] Building gene–peak graph...")
    links_df = pd.read_csv(args.links, sep="\t")
    peak_bed = scglue.genomics.read_bed(args.peaks)
    gene2peak = nx.Graph()
    for r in links_df.itertuples():
        gene2peak.add_edge(r.gene, r.peak_id, weight=getattr(r, 'score', 1.0))
    
    if args.gene2peak_out:
        print(f"[2.1] Saving gene–peak graph to: {args.gene2peak_out}")
        gene2peak_df = nx.to_pandas_edgelist(gene2peak)
        gene2peak_df.to_csv(args.gene2peak_out, index=False)
    genes = pd.Index(links_df['gene'].unique())
    peaks = pd.Index(links_df['peak_id'].unique())

    # 3. Build peak–TF graph via motif overlap
    print("[3/6] Building peak–TF graph...")
    motif_bed = scglue.genomics.read_bed(args.motif_bed)
    tfs = pd.Index(motif_bed["name"].unique())
    print(f"{len(tfs)} TFs retained for cis-ranking")
    peak2tf = scglue.genomics.window_graph(peak_bed, motif_bed, 0, right_sorted=True)
    peak2tf = peak2tf.edge_subgraph(e for e in peak2tf.edges if e[1] in tfs)
    if args.peak2tf_out:
        print(f"[3.1] Saving peak–TF graph to: {args.peak2tf_out}")
        peak2tf_df = nx.to_pandas_edgelist(peak2tf)
        peak2tf_df.to_csv(args.peak2tf_out, index=False)

    # 4. Rank TFs per gene with cis-regulatory evidence
    print("[4/6] Ranking TFs per gene...")
    region_lens = (links_df['chromEnd'] - links_df['chromStart']).values
    region_len = pd.Series(region_lens, index=links_df['peak_id'])
    gene2tf_rank = scglue.genomics.cis_regulatory_ranking(
        gene2peak, peak2tf,
        genes, peaks, tfs,
        region_lens=region_len,
        random_state=args.seed
    )
    gene2tf_rank.columns = gene2tf_rank.columns + "_glue"
    scglue.genomics.write_scenic_feather(gene2tf_rank, args.rankings)
    

    # Rank based on promoter                      
    adata = sc.read_loom(args.loom, sparse=True)
    scglue.data.get_gene_annotation(
        adata, gtf="gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz",
        gtf_by="gene_name"
    )
    gene_df = (
        adata.var.loc[:, ["chrom", "chromStart", "chromEnd", "strand"]]
        .dropna(subset=["chrom", "chromStart", "chromEnd"])
        .copy()
    )
    gene_df["chromStart"] = gene_df["chromStart"].astype(int) - 1   # 0-based BED
    gene_df["name"]       = gene_df.index                           # unique ID per row
    gene_df["score"]      = 0                                       # placeholder, keeps 6 cols 
    
    flank_bed = scglue.genomics.Bed(gene_df).strand_specific_start_site().expand(500, 500)
    flank2tf = scglue.genomics.window_graph(flank_bed, motif_bed, 0, right_sorted=True)
    if args.flank2tf_out:
        print(f"[3.1] Saving flank–TF graph to: {args.flank2tf_out}")
        flank2tf_df = nx.to_pandas_edgelist(flank2tf)
        flank2tf_df.to_csv(args.flank2tf_out, index=False)
    
    flank2tf_df = pd.read_csv(args.flank2tf_out)
    # Rebuild the graph from the DataFrame
    flank2tf = nx.from_pandas_edgelist(
        flank2tf_df,
        source="source",
        target="target"
    )
    
    genes   = gene_df.index.tolist()          # gene symbols
    regions = flank_bed.index.tolist()        # window IDs (they are NOT genes!)

    # gene2flank = nx.Graph([(g, g) for g in genes])
    gene2flank = nx.Graph((g, r) for g, r in zip(genes, regions))
    gene2tf_rank_supp = scglue.genomics.cis_regulatory_ranking(
        gene2flank, flank2tf, genes, regions, tfs,
        n_samples=0
    )
    gene2tf_rank_supp.columns = gene2tf_rank_supp.columns + "_supp"
    scglue.genomics.write_scenic_feather(gene2tf_rank_supp, args.rankings_prom)


    # 5. Create ctx_annotation.tsv for pyscenic ctx
    print("[5/6] Writing ctx annotation file...")
    pd.concat([
        pd.DataFrame({
            "#motif_id": tfs + "_glue",
            "gene_name": tfs
        }),
        pd.DataFrame({
            "#motif_id": tfs + "_supp",
            "gene_name": tfs
        })
    ]).assign(
        motif_similarity_qvalue=0.0,
        orthologous_identity=1.0,
        description="placeholder"
    ).to_csv(args.annotation, sep="\t", index=False)

    # 6. Prune draft network with pyscenic ctx
    cmd_ctx = [
        "pyscenic", "ctx",
        args.draft, args.rankings, args.rankings_prom,
        "--annotations_fname", args.annotation,
        "--expression_mtx_fname", args.loom,
        "--output", args.pruned,
        "--rank_threshold", "1000",
        "--min_genes", "1",
        "--num_workers", str(args.workers),
        "--cell_id_attribute", "CellID",
        "--gene_attribute", "Gene"
    ]
    print("[6/6] Running pyscenic ctx (pruning)...")
    subprocess.run(cmd_ctx, check=True)

    print("\nAll steps completed successfully.")


if __name__ == "__main__":
    main()
