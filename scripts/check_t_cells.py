import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("data/processed/pt011_clustered.h5ad")
raw = adata.raw.to_adata()

# Add the markers to the object for quick lookup
genes = ["CD3E", "CD3D", "CD8A", "CD4", "KLRF1", "FCGR3A", "CD14", "LYZ"]
present = [g for g in genes if g in raw.var_names]
print("Markers found:", present)

# Percent of cells expressing each marker, per cluster
expr = raw[:, present].to_df()
expr["cluster"] = adata.obs["cluster"].values
pct = expr.groupby("cluster").agg(lambda x: (x > 0).mean() * 100).round(0)
print("\nPercent of cells expressing each marker, by cluster:")
print(pct.to_string())

# T cells (CD3E or CD3D expressed) by group
t_cell = (expr["CD3E"] > 0) | (expr["CD3D"] > 0)
cd8 = expr["CD8A"] > 0
df = pd.DataFrame({
    "Cell_Source": adata.obs["Cell_Source"].values,
    "T_cell": t_cell.values,
    "CD8_T": (t_cell & cd8).values,
})
print("\nT cells and CD8 T cells per group:")
print(df.groupby("Cell_Source")[["T_cell", "CD8_T"]].sum().to_string())
