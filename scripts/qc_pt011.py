import numpy as np
import scanpy as sc

SRC = "data/raw/GSE266486_Patient_Aggregate_Dataset.h5ad"
DST = "data/processed/pt011.h5ad"

adata = sc.read_h5ad(SRC, backed="r")
sub = adata[adata.obs["Patient"] == "Pt011"].to_memory()
print("Pt011 cells:", sub.shape)

sample = sub.X.data[:5000] if hasattr(sub.X, "data") else np.asarray(sub.X).ravel()[:5000]
print("Integer counts:", bool(np.all(np.mod(sample, 1) == 0)))
print("Max value:", float(sample.max()))

sub.write_h5ad(DST)
print("Saved", DST)

sub.var["mt"] = sub.var_names.str.startswith("MT-")
print("Mitochondrial genes found:", int(sub.var["mt"].sum()))
sc.pp.calculate_qc_metrics(sub, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
print("\nOverall distribution:")
print(sub.obs[cols].describe(percentiles=[0.05, 0.5, 0.95]).round(2).to_string())

print("\nMedian by cell source:")
print(sub.obs.groupby("Cell_Source")[cols].median().round(2).to_string())

print("\nCells in each Cell_Source:")
print(sub.obs["Cell_Source"].value_counts().to_string())

# Count cells passing each suggested threshold (nothing is removed)
o = sub.obs
min_genes = o["n_genes_by_counts"] >= 500
max_genes = o["n_genes_by_counts"] <= 6000
low_mt = o["pct_counts_mt"] < 10
all_pass = min_genes & max_genes & low_mt

n = len(o)
print("\nCells passing each threshold:")
print(f"  >= 500 genes:       {int(min_genes.sum()):>6} ({min_genes.mean()*100:.1f}%)")
print(f"  <= 6000 genes:      {int(max_genes.sum()):>6} ({max_genes.mean()*100:.1f}%)")
print(f"  < 10% mito:         {int(low_mt.sum()):>6} ({low_mt.mean()*100:.1f}%)")
print(f"  all three together: {int(all_pass.sum()):>6} ({all_pass.mean()*100:.1f}%)")

sub.obs[["Cell_Source", "sample_id"] + cols].to_csv("docs/qc_metrics_pt011.csv")
print("\nPer-cell metrics saved to docs/qc_metrics_pt011.csv")
