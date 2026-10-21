import numpy as np
import scanpy as sc

SRC = "data/raw/GSE266486_Patient_Aggregate_Dataset.h5ad"
DST = "data/processed/pt006.h5ad"

# Load in backed mode, then pull only Pt006 into memory
adata = sc.read_h5ad(SRC, backed="r")
sub = adata[adata.obs["Patient"] == "Pt006"].to_memory()
print("Pt006 cells:", sub.shape)

# Check that the matrix holds raw integer counts, not normalized values
sample = sub.X.data[:5000] if hasattr(sub.X, "data") else np.asarray(sub.X).ravel()[:5000]
print("Integer counts:", bool(np.all(np.mod(sample, 1) == 0)))
print("Max value:", float(sample.max()))

sub.write_h5ad(DST)
print("Saved", DST)

# Mitochondrial genes (human symbols start with MT-)
sub.var["mt"] = sub.var_names.str.startswith("MT-")
print("Mitochondrial genes found:", int(sub.var["mt"].sum()))

sc.pp.calculate_qc_metrics(sub, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
print("\nOverall distribution:")
print(sub.obs[cols].describe(percentiles=[0.05, 0.5, 0.95]).round(2).to_string())

print("\nMedian by cell source:")
print(sub.obs.groupby("Cell_Source")[cols].median().round(2).to_string())

sub.obs[["Cell_Source", "sample_id"] + cols].to_csv("docs/qc_metrics_pt006.csv")
print("\nPer-cell metrics saved to docs/qc_metrics_pt006.csv")
