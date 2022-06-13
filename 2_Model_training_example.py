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
import itertools # This module implements a number of iterator building blocks
import networkx as nx
import pandas as pd # a fast, powerful, flexible and easy to use open source data analysis and manipulation tool
import scanpy as sc
import scglue
import seaborn as sns # a Python data visualization library based on matplotlib
from matplotlib import rcParams


# Parameters
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# Read preprocessed data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad")
atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad")
graph = nx.read_graphml("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/prior.graphml.gz")


# Configure data
# For each dataset to be integrated, we specify the probabilistic generative model to use
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)


# subset the prior graph to retain highly variable features only
graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))


# Build and train GLUE model
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, graph,
    fit_kws={"directory": "glue"}
)

glue.save("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/glue.dill")


# Check integration diagnostics
# To check whether the integration is reliable, we provide an “integration consistency score”, 
# which quantifies the consistency between the integration result and the guidance graph
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, graph,
    count_layers={"rna": "counts"}
)
dx

# _ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")


# Apply model for cell and feature embedding
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)


# Save the RNA and ATAC data
rna.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad", compression="gzip")
atac.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad", compression="gzip")


# Visualization
combined = anndata.concat([rna, atac])
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
# sc.tl.umap(combined)
# sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.65)


# Generate embeddings
feature_embeddings = glue.encode_graph(graph)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]
feature_embeddings.to_csv("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/feature_embeddings.csv", 
                            index=True, header=False)
