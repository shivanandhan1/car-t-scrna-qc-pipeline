# CAR-T single-cell QC and metadata validation

A quality-control and metadata-validation pipeline for public single-cell RNA-seq data from CAR-T cell therapy patients.

## Dataset

- **Source:** GEO accession GSE266486, from Kaminski J et al., *Blood* 2024 (PMID 38558106)
- **Patient analyzed:** Pt011 (12,244 cells: Bag, PBMC, CARNeg, CARPos)
- **Platform:** Illumina NovaSeq 6000 (GPL24676)
- **Library type:** [10x Chromium, confirmed via SRA Run Selector]

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
- CARPos 88.7% kept, CARNeg 96.9%, PBMC 97.9%, Bag 56.7%
- Metadata validation: 15 of 15 checks passed

## Trade-offs

- **The mitochondrial filter removes 43% of the Bag group.** The Bag group has a median mitochondrial rate of 10.4%, nearly double the other groups. This may reflect biology, such as stress during processing, rather than poor data quality. The filter is applied uniformly, so it changes the group mix. Any comparison across groups should account for that.
- **Validation rules are project-defined.** They are not full CELLxGENE schema compliance.
- **One alias was applied.** The source labels the Bag group `Pt011_Bag` in `sample_id` and `Pt011_Product` in `Formatted_Sample_Name`. I treated them as one group. [Confirmed from paper methods / still to confirm]

## Not in this project

- No raw sequencing reads (the authors submitted them to dbGaP, and the raw tar was not used)
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
