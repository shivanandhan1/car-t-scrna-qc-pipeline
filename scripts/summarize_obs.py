import scanpy as sc

adata = sc.read_h5ad("data/raw/GSE266486_Patient_Aggregate_Dataset.h5ad", backed="r")
obs = adata.obs

print("\nSample IDs and cell counts:")
print(obs["sample_id"].value_counts().to_string())

print("\nPatient IDs:")
print(obs["Patient"].value_counts().to_string())

print("\nCell source:")
print(obs["Cell_Source"].value_counts().to_string())

print("\nMissing values per column:")
print(obs.isna().sum().to_string())
