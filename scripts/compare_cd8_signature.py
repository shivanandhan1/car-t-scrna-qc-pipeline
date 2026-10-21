import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import mannwhitneyu

adata = sc.read_h5ad("data/processed/pt011_clustered.h5ad")

# CD8 T-cell clusters from the marker check, comparing only the two CAR groups
cd8_clusters = ["3", "4", "8", "9", "13"]
sub = adata[adata.obs["cluster"].isin(cd8_clusters)
            & adata.obs["Cell_Source"].isin(["CARPos", "CARNeg"])].copy()

# Use the full log-normalized data, not the 2,000 variable genes
raw = sub.raw.to_adata()

up_genes = ["KIR3DL2", "CD160", "KLRD1", "CCL5", "XCL1", "CCR9"]
down_genes = ["SELL", "CD28"]
for name, genes in [("up", up_genes), ("down", down_genes)]:
    missing = [g for g in genes if g not in raw.var_names]
    print(f"{name} genes missing:", missing if missing else "none")

sc.tl.score_genes(raw, [g for g in up_genes if g in raw.var_names], score_name="bystander_up")
sc.tl.score_genes(raw, [g for g in down_genes if g in raw.var_names], score_name="naive_down")
raw.obs["net_score"] = raw.obs["bystander_up"] - raw.obs["naive_down"]

df = raw.obs[["cluster", "Cell_Source", "bystander_up", "naive_down", "net_score"]]

print("\nCell counts and mean scores by cluster and group:")
print(df.groupby(["cluster", "Cell_Source"], observed=True).agg(
    n=("net_score", "size"),
    bystander_up=("bystander_up", "mean"),
    naive_down=("naive_down", "mean"),
    net_score=("net_score", "mean"),
).round(3).to_string())

# Compare CARNeg with CARPos within each cluster
print("\nCARNeg minus CARPos, per cluster:")
rows = []
for c, g in df.groupby("cluster", observed=True):
    a = g.loc[g["Cell_Source"] == "CARNeg", "net_score"]
    b = g.loc[g["Cell_Source"] == "CARPos", "net_score"]
    p = mannwhitneyu(a, b).pvalue if len(a) >= 20 and len(b) >= 20 else np.nan
    rows.append((c, len(a), len(b), round(a.mean() - b.mean(), 3), p))
print(pd.DataFrame(rows, columns=["cluster", "n_CARNeg", "n_CARPos", "diff", "p_value"]).to_string(index=False))

# Pooled comparison across all CD8 clusters
a = df.loc[df["Cell_Source"] == "CARNeg", "net_score"]
b = df.loc[df["Cell_Source"] == "CARPos", "net_score"]
print(f"\nPooled: CARNeg {a.mean():.3f} vs CARPos {b.mean():.3f}, "
      f"p = {mannwhitneyu(a, b).pvalue:.3g}")
