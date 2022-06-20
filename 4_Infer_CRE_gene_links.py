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


# %%
import functools
import gc
import itertools
import operator
import os
import sys
from math import ceil

import anndata
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import scipy.stats
import seaborn as sns
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

import scglue
import utils


# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


np.random.seed(0)


# Get arguments
out_dir = sys.argv[1] # output directory
if not os.path.exists(out_dir):
    print('Creating the directory ' + out_dir + ' ...\n')
    os.makedirs(out_dir)


# Load data
rna = anndata.read_h5ad(f"{out_dir}/rna_preprocessed.h5ad")
atac = anndata.read_h5ad(f"{out_dir}/atac_preprocessed.h5ad")


# %%
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names).query("dcq_highly_variable"))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names).query("dcq_highly_variable"))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)
flanks = tss.expand(500, 500)


# Load graph
dist_graph = nx.read_graphml(f"{out_dir}/dist.graphml.gz")
# Serves as genomic windows

dist = biadjacency_matrix(dist_graph, genes.index, peaks.index, weight = "dist", dtype = np.float32)

# # %%
# corr = biadjacency_matrix(
#     scglue.data.metacell_corr(
#         *rna, *atac, "X_pca", n_meta=200,
#         skeleton=dist_graph.subgraph([*genes.index, *peaks.index]), method="spr"
#     ), genes.index, peaks.index, weight="corr", dtype=np.float32
# )


# Load feature embeddings
feature_embeddings = pd.read_csv(f"{out_dir}/feature_embeddings.csv", 
                                  header = None, index_col = 0)


# %%
glue_graph = scglue.genomics.regulatory_inference(
    feature_embeddings.index,
    feature_embeddings.to_numpy(),
    dist_graph.subgraph([*genes.index, *peaks.index]),
    alternative="greater", random_state=0
)
glue = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="score", dtype=np.float32)
qval = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="qval", dtype=np.float32)


# %%
window = biadjacency_matrix(
    dist_graph, genes.index, peaks.index, weight=None, dtype=np.float32
).tocoo()
dist = window.multiply(dist.toarray())
# corr = window.multiply(corr.toarray())
glue = window.multiply(glue.toarray())
qval = window.multiply(qval.toarray())



for mat in (dist, glue, qval):
    assert np.all(window.row == mat.row)
    assert np.all(window.col == mat.col)


# %%
gene_peak_conn = pd.DataFrame({
    "gene": genes.index[window.row],
    "peak": peaks.index[window.col],
    "dist": dist.data.astype(int), # dist : two choices
    # "corr": corr.data,
    "glue": glue.data, # glue : two choices
    "qval": qval.data # qval ; two choices
})
gene_peak_conn.to_pickle(f"{out_dir}/gene_peak_conn.pkl.gz")
gene_peak_df = pd.DataFrame(gene_peak_conn)
gene_peak_df.to_csv(f"{out_dir}/gene_peak_conn.csv", 
                            index = False, header = True)


# # Filtering gene-peak connection
# # %%
# qval_cutoff = 0.05
# _ = sns.scatterplot(
#     x="glue", y="qval", data=gene_peak_conn.sample(n=2000),
#     edgecolor=None, s=5, alpha=0.1
# ).axhline(y=qval_cutoff, c="darkred", ls="--")
# 
# 
# # %%
# gene_peak_conn_glue = gene_peak_conn.query(f"qval < {qval_cutoff}")
# gene_peak_conn_glue.shape[0]
# 
# 
# # %%
# frac_pos = gene_peak_conn_glue.shape[0] / gene_peak_conn.shape[0]
# frac_pos
# 
# # %%
# glue_cutoff = gene_peak_conn.query(f"qval < {qval_cutoff}")["glue"].min()
# glue_cutoff
