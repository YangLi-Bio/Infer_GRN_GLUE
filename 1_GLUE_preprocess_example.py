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


# Import modules
import anndata
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams


# Parameters
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)


# Load data
rna = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Chen-2019-RNA.h5ad")
# rna

atac = anndata.read_h5ad("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/Chen-2019-ATAC.h5ad")
# atac


# Preprocess scRNA-seq data
rna.X, rna.X.data
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.umap(rna, color="cell_type")


# Preprocess scATAC-seq data
atac.X, atac.X.data
scglue.data.lsi(atac, n_components=100, n_iter=15)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)
sc.pl.umap(atac, color="cell_type")


# Obtain genomic coordinates
rna.var.head()
scglue.data.get_gene_annotation(
    rna, gtf="/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Annotations/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
atac.var_names[:5]
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1])
atac.var["chromEnd"] = split.map(lambda x: x[2])
atac.var.head()


# Graph construction
graph = scglue.genomics.rna_anchored_prior_graph(rna, atac)
graph
graph.number_of_nodes(), graph.number_of_edges()


# Graph node covers all omic features
all(graph.has_node(gene) for gene in rna.var_names), \
all(graph.has_node(peak) for peak in atac.var_names)


# Edge attributes contain weights and signs
for _, e in zip(range(5), graph.edges):
    print(f"{e}: {graph.edges[e]}")


# Each node has a self-loop
all(graph.has_edge(gene, gene) for gene in rna.var_names), \
all(graph.has_edge(peak, peak) for peak in atac.var_names)


# Graph is symmetric
all(graph.has_edge(j, i) for i, j, _ in graph.edges)
atac.var.head()


rna.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/rna_preprocessed.h5ad", compression="gzip")
atac.write("/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/atac_preprocessed.h5ad", compression="gzip")
nx.write_graphml(graph, "/fs/ess/PCON0022/liyang/STREAM/benchmarking/GLUE/Example/prior.graphml.gz")
