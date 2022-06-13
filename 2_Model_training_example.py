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
import os
import re

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams

import scglue


# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# %%
PRIOR = os.environ.get("PRIOR", "d")
SEED = int(os.environ.get("SEED", "0"))


# Load data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad")
rna.X, rna.obs, rna.var
atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad")
atac.X, atac.obs, atac.var


# %%
scglue.models.configure_dataset(rna, "NB", use_highly_variable=True, use_rep="X_pca")
scglue.models.configure_dataset(atac, "NB", use_highly_variable=True, use_rep="X_lsi")
# %% tags=[]
glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(prior.nodes),
    random_seed=SEED
)
glue.compile()
glue.fit(
    {"rna": rna, "atac": atac},
    prior, edge_weight="weight", edge_sign="sign",
    align_burnin=np.inf, safe_burnin=False,
    directory="/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/"
)
glue.save("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Pretrain.dill")


# %%
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
scglue.data.estimate_balancing_weight(
    rna, atac, use_rep="X_glue"
)

# %%
scglue.models.configure_dataset(rna, "NB", use_highly_variable=True, use_rep="X_pca", use_dsc_weight="balancing_weight")
scglue.models.configure_dataset(atac, "NB", use_highly_variable=True, use_rep="X_lsi", use_dsc_weight="balancing_weight")

# %%
glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(prior.nodes),
    random_seed=SEED
)
glue.adopt_pretrained_model(scglue.models.load_model(
    "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Pretrain.dill"
))
glue.compile()
glue.fit(
    {"rna": rna, "atac": atac},
    prior, edge_weight="weight", edge_sign="sign",
    directory="/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/"
)
glue.save("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Final.dill")


# %%
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

# %%
combined = anndata.AnnData(
    obs=pd.concat([rna.obs, atac.obs], join="inner"),
    obsm={"X_glue": np.concatenate([rna.obsm["X_glue"], atac.obsm["X_glue"]])}
)

# %%
sc.pp.neighbors(combined, n_pcs=50, use_rep="X_glue", metric="cosine")
# sc.tl.umap(combined)


# %%
rna.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_glue.h5ad", compression="gzip")
atac.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_glue.h5ad", compression="gzip")
combined.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/combined_glue.h5ad", compression="gzip")


# %% [markdown]
# ## Feature embeddings

# %%
feature_embeddings = pd.DataFrame(
    glue.encode_graph(prior, "weight", "sign"),
    index=glue.vertices
)
feature_embeddings.iloc[:5, :5]

# %%
feature_embeddings.to_csv("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/feature_embeddings.csv", 
                            index=True, header=False)
