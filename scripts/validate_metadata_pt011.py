import re
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("data/processed/pt011_filtered.h5ad")
obs = adata.obs.copy()
for col in obs.columns:
    if isinstance(obs[col].dtype, pd.CategoricalDtype):
        obs[col] = obs[col].astype(str)
results = []

def check(name, passed, detail=""):
    results.append((name, "PASS" if passed else "FAIL", detail))

# 1. Required fields are present and complete
required = ["sample_id", "donor", "origin", "Patient", "Cell_Source",
            "Formatted_Sample_Name", "gem_well"]
for col in required:
    if col not in obs.columns:
        check(f"{col} present", False, "column missing")
    else:
        missing = int(obs[col].isna().sum())
        check(f"{col} present and complete", missing == 0, f"{missing} missing")

# 2. Cell_Source uses only the expected labels
allowed_sources = {"Bag", "PBMC", "CARPos", "CARNeg"}
unexpected = sorted(set(obs["Cell_Source"]) - allowed_sources)
check("Cell_Source uses expected labels", not unexpected,
      f"unexpected: {unexpected}" if unexpected else "all expected")

# 3. Only the target patient remains
patients = sorted(obs["Patient"].unique())
check("Single patient (Pt011)", patients == ["Pt011"], f"found: {patients}")

# 4. Donor matches Patient
mismatch = int((obs["donor"] != obs["Patient"]).sum())
check("donor matches Patient", mismatch == 0, f"{mismatch} mismatches")

# 5. sample_id follows Patient_Cell_Source
expected_id = obs["Patient"] + "_" + obs["Cell_Source"]
mismatch = int((obs["sample_id"] != expected_id).sum())
check("sample_id = Patient_Cell_Source", mismatch == 0, f"{mismatch} mismatches")

# 6. Formatted_Sample_Name matches sample_id, using a documented alias
ALIASES = {"Pt011_Bag": "Pt011_Product"}  # source labels Bag as Product; confirm in paper methods
expected_name = obs["sample_id"].replace(ALIASES)
mismatch = int((obs["Formatted_Sample_Name"] != expected_name).sum())
alias_count = int(obs["sample_id"].isin(ALIASES.keys()).sum())
check("Formatted_Sample_Name matches sample_id (alias applied)", mismatch == 0,
      f"{mismatch} mismatches after alias; alias used for {alias_count} cells")

# 7. Cell barcodes are unique and well-formed
dupes = int(obs.index.duplicated().sum())
check("cell barcodes unique", dupes == 0, f"{dupes} duplicates")
bad_format = int((~obs.index.to_series().str.match(r"^[ACGT]+-\d+$")).sum())
check("cell barcode format (e.g. AAACCTGAGGAATTAC-1)", bad_format == 0,
      f"{bad_format} malformed")

# 8. gem_well is a positive integer
gw = pd.to_numeric(obs["gem_well"], errors="coerce")
bad_gw = int((gw.isna() | (gw < 1)).sum())
check("gem_well is a positive integer", bad_gw == 0, f"{bad_gw} invalid")

# Report
report = pd.DataFrame(results, columns=["check", "result", "detail"])
print(report.to_string(index=False))
print(f"\n{(report['result'] == 'PASS').sum()} of {len(report)} checks passed")

with open("docs/metadata_validation_pt011.md", "w") as f:
    f.write("# Metadata validation: Pt011 filtered data\n\n")
    f.write("| Check | Result | Detail |\n|---|---|---|\n")
    for _, row in report.iterrows():
        f.write(f"| {row['check']} | {row['result']} | {row['detail']} |\n")
print("\nReport saved to docs/metadata_validation_pt011.md")
