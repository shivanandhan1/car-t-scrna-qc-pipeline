import scanpy as sc

MIN_GENES, MAX_GENES, MAX_MT_PCT = 500, 6000, 15

adata = sc.read_h5ad("data/processed/pt011.h5ad")
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

keep = (
    (adata.obs["n_genes_by_counts"] >= MIN_GENES)
    & (adata.obs["n_genes_by_counts"] <= MAX_GENES)
    & (adata.obs["pct_counts_mt"] < MAX_MT_PCT)
)

filtered = adata[keep].copy()
filtered.write_h5ad("data/processed/pt011_filtered.h5ad")

print(f"Before: {adata.n_obs} cells")
print(f"After:  {filtered.n_obs} cells ({filtered.n_obs / adata.n_obs * 100:.1f}% kept)")
print("\nKept per cell source:")
before = adata.obs["Cell_Source"].value_counts()
after = filtered.obs["Cell_Source"].value_counts()
table = [(src, after.get(src, 0), before[src]) for src in before.index]
for src, kept, total in table:
    print(f"  {src:<8} {kept:>6} of {total:<6} ({kept / total * 100:.1f}%)")
print("\nSaved data/processed/pt011_filtered.h5ad")
