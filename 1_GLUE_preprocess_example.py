# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---


# Modules
import functools #  higher-order functions that work on other functions
import itertools #  implements a number of iterator building blocks inspired by constructs 
# from APL, Haskell, and SML

import os # The design of all built-in operating system dependent modules

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix

import scglue


# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# Load data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Chen-2019-RNA.h5ad")
rna.X, rna.obs, rna.var
atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Chen-2019-ATAC.h5ad")
atac.X, atac.obs, atac.var


# # Set cell types and chromosomes
# rna.obs["cell_type"].cat.categories
# 
# # %%
# used_cts = {
#     "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TEM_1", "CD8 TEM_2",
#     "CD14 Mono", "CD16 Mono", "Memory B", "Naive B"
# }  # To match cell types covered in PC Hi-C
# used_chroms = {f"chr{x}" for x in range(1, 23)}.union({"chrX"})
# 
# # %%
# rna = rna[
#     [item in used_cts for item in rna.obs["cell_type"]],
#     [item in used_chroms for item in rna.var["chrom"]]
# ]

sc.pp.filter_genes(rna, min_counts=1)
rna.obs_names += "-RNA"

# rna
# 
# # %%
# atac = atac[
#     [item in used_cts for item in atac.obs["cell_type"]],
#     [item in used_chroms for item in atac.var["chrom"]]
# ]

sc.pp.filter_genes(atac, min_counts=1)
atac.obs_names += "-ATAC"

# atac


# Genes, CREs, promoters, and TSSs
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)


# Preprocess RNA
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=6000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, n_comps=100, use_highly_variable=True, svd_solver="auto")
sc.pp.neighbors(rna, n_pcs=100, metric="cosine")
# sc.tl.umap(rna)

# %%
rna.X = rna.layers["counts"]
del rna.layers["counts"]


# Preprocess ATAC
scglue.data.lsi(atac, n_components=100, use_highly_variable=False, n_iter=15)
sc.pp.neighbors(atac, n_pcs=100, use_rep="X_lsi", metric="cosine")
# sc.tl.umap(atac)


# Generate prior overlap graph
overlap_graph = scglue.genomics.window_graph(
    genes.expand(2000, 0), peaks, 0,
    attr_fn=lambda l, r, d: {
        "weight": 1.0,
        "type": "overlap"
    }
)
overlap_graph = nx.DiGraph(overlap_graph)
overlap_graph.number_of_edges()


# Build distance graph
dist_graph = scglue.genomics.window_graph(
    promoters, peaks, 150000,
    attr_fn=lambda l, r, d: {
        "dist": abs(d),
        "weight": scglue.genomics.dist_power_decay(abs(d)),
        "type": "dist"
    }
)
dist_graph = nx.DiGraph(dist_graph)
dist_graph.number_of_edges()


# %%
rna.var["dcq_highly_variable"] = rna.var["highly_variable"]
rna.var["dcq_highly_variable"].sum()

# %%
hvg_reachable = scglue.graph.reachable_vertices(dist_graph, rna.var.query("dcq_highly_variable").index)

# %%
atac.var["dcq_highly_variable"] = [item in hvg_reachable for item in atac.var_names]
atac.var["dcq_highly_variable"].sum()


# Save the preprocessed data
rna.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad", compression="gzip")
atac.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad", compression="gzip")
nx.write_graphml(overlap_graph, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/overlap.graphml.gz")
nx.write_graphml(dist_graph, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/dist.graphml.gz")
