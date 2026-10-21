import pandas as pd

df = pd.read_csv("docs/qc_metrics_pt011.csv", index_col=0)

for mt in [10, 15, 20]:
    keep = (
        (df["n_genes_by_counts"] >= 500)
        & (df["n_genes_by_counts"] <= 6000)
        & (df["pct_counts_mt"] < mt)
    )
    print(f"\nmito < {mt}%  (genes 500-6000 applied)")
    summary = df.assign(keep=keep).groupby("Cell_Source")["keep"].agg(["sum", "count", "mean"])
    summary["mean"] = (summary["mean"] * 100).round(1)
    summary.columns = ["kept", "total", "percent_kept"]
    print(summary.to_string())
    print(f"  overall kept: {int(keep.sum())} of {len(df)}")
