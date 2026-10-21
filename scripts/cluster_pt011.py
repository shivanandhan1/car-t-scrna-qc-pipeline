import matplotlib
matplotlib.use("Agg")
import pandas as pd
import scanpy as sc

sc.settings.figdir = "docs"

adata = sc.read_h5ad("data/processed/pt011_filtered.h5ad")
adata.layers["counts"] = adata.X.copy()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Keep the full log data for marker plots, then reduce to variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()

# Dimensionality reduction and clustering
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="cluster")

print("Clusters:", adata.obs["cluster"].nunique())
print("\nCells per cluster:")
print(adata.obs["cluster"].value_counts().sort_index().to_string())

print("\nCluster by Cell_Source:")
print(pd.crosstab(adata.obs["cluster"], adata.obs["Cell_Source"]).to_string())

# Marker genes for cell-type assignment
markers = {
    "CD8 T": ["CD8A", "CD8B"],
    "CD4 T": ["CD4", "IL7R"],
    "Cytotoxic / NK": ["GZMB", "NKG7", "KLRD1"],
    "Naive T": ["SELL", "CCR7", "TCF7"],
    "B cell": ["MS4A1", "CD79A"],
    "Monocyte": ["CD14", "LYZ"],
    "Proliferating": ["MKI67", "TOP2A"],
}
available = set(adata.raw.var_names)
marker_list = [g for genes in markers.values() for g in genes if g in available]
missing = [g for genes in markers.values() for g in genes if g not in available]
print("\nMarkers not found in data:", missing if missing else "none")

sc.pl.umap(adata, color=["cluster", "Cell_Source"], save="_clusters.png", show=False)
sc.pl.dotplot(adata, marker_list, groupby="cluster", use_raw=True,
              save="_markers.png", show=False)

adata.write_h5ad("data/processed/pt011_clustered.h5ad")
print("\nSaved data/processed/pt011_clustered.h5ad")
print("Figures saved to docs/figures/")
