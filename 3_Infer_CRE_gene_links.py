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
import functools # higher-order functions that work on other functions
import gc # provides an interface to the optional garbage collector
import itertools #  implements a number of iterator building blocks inspired by constructs from APL, Haskell, and SML
import operator # divides the operators in the following groups: Arithmetic operators; 
# Assignment operators; Comparison operators; Logical operators ...

import os
from math import ceil

import anndata
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd # a fast, powerful, flexible and easy to use open source data analysis and manipulation tool
import scanpy as sc
import scipy.sparse # Sparse matrices
import scipy.stats # contains a large number of probability distributions, summary and frequency statistics, 
# correlation functions and statistical tests

import seaborn as sns # a Python data visualization library based on matplotlib
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform # Distance matrix computation from a collection 
# of raw observation vectors stored in a rectangular array

import scglue
import utils # a collection of small Python functions and classes which make common patterns shorter and easier.


# %%
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


np.random.seed(0)


# Load data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad")
atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad")


# %%
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names).query("dcq_highly_variable"))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names).query("dcq_highly_variable"))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)
flanks = tss.expand(500, 500)
