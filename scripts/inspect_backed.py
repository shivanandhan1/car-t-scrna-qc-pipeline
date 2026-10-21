import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE266486_Patient_Aggregate_Dataset.h5ad", backed="r")
print(adata.shape)
print(adata.obs.columns.tolist())
print(adata.obs.iloc[:5])
