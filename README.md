# CAR-T single-cell QC and metadata validation

A quality-control and metadata-validation pipeline for public single-cell RNA-seq data from CAR-T cell therapy patients.

## Dataset

- **Source:** GEO accession GSE266486, from Kaminski J et al., *Blood* 2024 (PMID 38558106)
- **Patient analyzed:** Pt011 (12,244 cells: CAR T infusion product, PBMC, CARNeg, CARPos)
- **Platform:** Illumina NovaSeq 6000 (GPL24676)
- **Library:** 10x Genomics Next GEM 5′ v2, with scRNA-seq and scTCR-seq (per GEO sample GSM8248030)
- **Processing by authors:** Cell Ranger multi v7.0.0, GRCh38 reference

The full dataset covers six patients. This project analyzes one, because it has the CAR-positive and CAR-negative groups needed for the paper's central comparison and fits on an 8 GB laptop.

## Pipeline

1. **Subset** to one patient and check that counts are raw integers.
2. **QC metrics:** genes per cell, UMIs per cell, and mitochondrial read percentage.
3. **Threshold selection:** tested mitochondrial cutoffs of 10%, 15%, and 20% and measured the effect on each cell group.
4. **Filter:** 500 to 6,000 genes per cell and less than 15% mitochondrial reads.
5. **Metadata validation:** 15 project-defined checks on identity, labels, and barcodes.

Scripts are in `scripts/`. Every step writes a record to `docs/` or `LINEAGE.md`.

## Results

- 12,244 cells in, 8,834 kept (72.1%)
- CARPos 88.7% kept, CARNeg 96.9%, PBMC 97.9%, infusion product (Bag) 56.7%
- Metadata validation: 15 of 15 checks passed

## Trade-offs

- **The mitochondrial filter removes 43% of the infusion product cells.** The Bag group has a median mitochondrial rate of 10.4%, about twice the other groups. This may reflect biology, such as stress during processing, rather than poor data quality. The filter is applied uniformly, so it changes the group mix. Any analysis of the infusion product itself is based on the cells that pass this filter.
- **Validation rules are project-defined.** They are not full CELLxGENE schema compliance.
- **One label alias was applied.** The source labels the infusion product as `Pt011_Bag` in `sample_id` and `Pt011_Product` in `Formatted_Sample_Name`. Confirmed as the same group in GEO sample GSM8248030.
- **Library chemistry is 5′ v2.** It captures TCR data alongside gene expression. The TCR data is missing for 46% of cells, so any TCR analysis covers only part of the dataset.

## Not in this project

- No raw sequencing reads. The authors state these are being submitted to dbGaP, and a search of SRA returned no items on 2026-10-21. The raw tar was not used.
- No clustering, cell-type annotation, or differential expression. These are the natural next steps.
- No warehouse or orchestration layer. The pipeline is in Python and scanpy, not dbt or Airflow.

## Reproducing

```bash
python -m venv .venv && source .venv/bin/activate
pip install scanpy anndata pandas duckdb
python scripts/qc_pt011.py
python scripts/pass_rates_pt011.py
python scripts/filter_pt011.py
python scripts/validate_metadata_pt011.py
```

Raw data lives in `data/raw/` and is not committed. See `LINEAGE.md` for the download source and checksum.
