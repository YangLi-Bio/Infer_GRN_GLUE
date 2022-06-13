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
import anndata
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams


# Parameters
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# Read preprocessed data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad")
atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad")
graph = nx.read_graphml("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/prior.graphml.gz")


# Configure data
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))


# Build and train GLUE model
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "glue"}
)

glue.save("glue.dill")
# glue = scglue.models.load_model("glue.dill")


# Check integration diagnostics
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, graph,
    count_layers={"rna": "counts"}
)
dx

_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")


# Apply model for cell and feature embedding
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = anndata.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.65)

feature_embeddings = glue.encode_graph(graph)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]
feature_embeddings.to_csv("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/feature_embeddings.csv", index=True, header=False)


# %%
corr = biadjacency_matrix(
    utils.metacell_corr(
        rna, atac, "X_pca", n_meta=200,
        skeleton=dist_graph.subgraph([*genes.index, *peaks.index]), method="spr"
    ), genes.index, peaks.index, weight="corr", dtype=np.float32
)


# Obtain genes, CREs, and promoters
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names).query("dcq_highly_variable"))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names).query("dcq_highly_variable"))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)
flanks = tss.expand(500, 500)


# Load modules
import functools
import itertools
import os
import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
import scglue


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
nx.write_graphml(dist_graph, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/dist.graphml.gz")


# %%
glue_graph = scglue.genomics.regulatory_inference(
    feature_embeddings.index, # Check this
    [feature_embedding.to_numpy() for feature_embedding in feature_embeddings],
    dist_graph.subgraph([*genes.index, *peaks.index]),
    alternative="greater", random_state=0
)


glue = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="score", dtype=np.float32)
qval = biadjacency_matrix(glue_graph, genes.index, peaks.index, weight="qval", dtype=np.float32)


dist = window.multiply(dist.toarray())
corr = window.multiply(corr.toarray())
glue = window.multiply(glue.toarray())
qval = window.multiply(qval.toarray())


# %%
gene_peak_conn = pd.DataFrame({
    "gene": genes.index[window.row],
    "peak": peaks.index[window.col],
    "dist": dist.data.astype(int),
    "corr": corr.data,
    "glue": glue.data,
    "qval": qval.data
})


gene_peak_conn.to_pickle("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/gene_peak_conn.pkl.gz")
