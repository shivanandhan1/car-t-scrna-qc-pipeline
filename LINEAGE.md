# Data lineage

Record of every dataset, file, and transformation in this project, in order.

## Source dataset

- **Accession:** GSE266486 (GEO)
- **Title:** B cell directed CAR-T cell therapy results in activation of CD8+ cytotoxic CARneg bystander T cells in both non-human primates and patients (human)
- **Source URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266486
- **Publication:** Kaminski J, Fleming RA, et al. *Blood* 2024;144(1):46-60. PMID 38558106. DOI 10.1182/blood.2023022717
- **Platform:** GPL24676 (Illumina NovaSeq 6000)
- **Submitted:** 2024-05-02. **Last GEO update:** 2024-08-13
- **Raw data status:** The authors state raw data is being submitted to dbGaP. As of 2026-10-21, a search of SRA for BioProject PRJNA1107270 returned no items. Raw reads were not used in this project.

## Verified details

| Item | Value | Source |
|---|---|---|
| Library chemistry | 10x Genomics Next GEM 5′ v2, scRNA-seq and scTCR-seq | GEO sample GSM8248030, Extraction protocol field |
| Processing | Cell Ranger multi v7.0.0, reference GRCh38 (Ensembl v107) | GEO sample GSM8248030, Data processing field |
| Bag group | CAR T infusion product. Source name "CAR T Product"; alt sample name Pt011_Bag | GEO sample GSM8248030 (titled Pt011_Product) |
| CAR construct | Anti-CD19 CAR, pediatric patients with B-ALL | GEO sample GSM8248030, Extraction protocol field |
| Cell barcodes | 16-base 10x format with -1 suffix, e.g. AAACCTGAGGAATTAC-1 | Observed in the .h5ad |

## Downloaded files

| File | Size (GEO) | Downloaded | Used | SHA-256 |
|---|---|---|---|---|
| GSE266486_Patient_Aggregate_Dataset.h5ad | 1.5 Gb | 2026-10-21 | Yes | `0a42f45e536c9ea8ac6b6e0f9f28145c9b426c9ba9b1bf85dee396629378cb3e` |
| GSE266486_RAW.tar | 1.0 Gb | Not downloaded | No | n/a |
| GSE266486_Patient_filtered_contig_annotations.csv.gz | 9.5 Mb | Not downloaded | No | n/a |

Stored in `data/raw/`. Raw files are never edited.

## Dataset contents (as read from the .h5ad)

- **Full object:** 82,416 cells x 20,669 genes
- **Patients:** Pt006 through Pt011, with 21 sample IDs matching the 21 patient samples on the GEO page. No non-human primate samples appear in this file.
- **Cells per patient:** Pt006 11,075; Pt007 29,437; Pt008 18,220; Pt009 8,795; Pt010 2,645; Pt011 12,244
- **Cell_Source labels:** Bag (CAR T infusion product), PBMC, CARPos, CARNeg, TCells
- **Missing metadata:** None in the fields used for validation. `raw_clonotype_id` is empty for 38,080 cells (46%).

## Patient selection

- **Pt006 (explored, not used):** 11,075 cells with CARPos, PBMC, and TCells. No CARNeg group, so the paper's central CARNeg comparison cannot be done with this patient.
- **Pt011 (selected):** 12,244 cells with Bag, PBMC, CARNeg, and CARPos. Covers the CAR-positive versus CAR-negative comparison and fits in 8 GB of RAM.

## Processing steps

### 1. Subset to Pt011
- **Script:** `scripts/qc_pt011.py`
- **Input:** `data/raw/GSE266486_Patient_Aggregate_Dataset.h5ad`
- **Output:** `data/processed/pt011.h5ad` (12,244 cells)
- **Checks:** Matrix holds integer counts (True); max value 219; 13 mitochondrial genes (`MT-` prefix)

### 2. QC metrics (no cells removed)
- **Metrics saved to:** `docs/qc_metrics_pt011.csv`
- **Pt011 medians:** 2,129.5 genes per cell; 4,825.5 UMIs per cell; 5.06% mitochondrial reads
- **Median mitochondrial % by source:** Bag 10.42; CARNeg 2.56; CARPos 4.44; PBMC 2.29

### 3. Threshold selection
- **Gene range:** 500 to 6,000 genes per cell. Pass rate 99.9% at the upper limit, 89.1% at the lower.
- **Mitochondrial cutoff tested:**

| Cutoff | Bag kept | CARNeg kept | CARPos kept | PBMC kept | Overall kept |
|---|---|---|---|---|---|
| 10% | 48.0% | 94.3% | 79.7% | 95.7% | 8,030 of 12,244 |
| **15% (chosen)** | **56.7%** | **96.9%** | **88.7%** | **97.9%** | **8,834 of 12,244** |
| 20% | 63.8% | 97.9% | 91.9% | 98.3% | 9,417 of 12,244 |

- **Script:** `scripts/pass_rates_pt011.py`
- **Decision:** 15% mitochondrial cutoff, chosen to keep most CARPos cells while removing high-mitochondrial cells.
- **Known trade-off:** The Bag group (CAR T infusion product) loses 43% of its cells at this cutoff. Its higher mitochondrial reads may reflect real biology from processing or stress. Any analysis of the product itself is based on the cells that pass this filter. Flag this in any write-up.

### 4. Filtering
- **Script:** `scripts/filter_pt011.py`
- **Filters:** 500 to 6,000 genes per cell; less than 15% mitochondrial reads
- **Output:** `data/processed/pt011_filtered.h5ad`
- **Result:** 8,834 of 12,244 cells kept (72.1%)

| Cell_Source | Kept | Total | Percent kept |
|---|---|---|---|
| Bag | 4,228 | 7,457 | 56.7% |
| PBMC | 2,129 | 2,174 | 97.9% |
| CARNeg | 1,886 | 1,947 | 96.9% |
| CARPos | 591 | 666 | 88.7% |

### 5. Metadata validation
- **Script:** `scripts/validate_metadata_pt011.py`
- **Output:** `docs/metadata_validation_pt011.md`
- **Result:** 15 of 15 checks passed
- **Checks:** required fields complete; Cell_Source labels expected; single patient (Pt011); donor matches Patient; sample_id equals Patient_Cell_Source; Formatted_Sample_Name matches sample_id; cell barcodes unique and well-formed; gem_well a positive integer
- **Alias applied:** `Pt011_Bag` in `sample_id` is labeled `Pt011_Product` in `Formatted_Sample_Name` for 4,228 cells. Confirmed as the same group: GEO sample GSM8248030 is titled Pt011_Product and has alt sample name Pt011_Bag.
- **Scope:** These are project-defined rules, not full CELLxGENE schema compliance.

### 6. Clustering, annotation, and bystander check (exploratory)
- Scripts: `scripts/cluster_pt011.py`, `scripts/check_t_cells.py`, `scripts/compare_cd8_signature.py`
- Output: `data/processed/pt011_clustered.h5ad`, figures in `docs/figures/`
- Clustering: 2,000 highly variable genes, 30 PCs, Leiden resolution 0.5, 15 clusters
- Annotation: CD8 T (3, 4, 8, 9, 13); CD4 T (0, 1, 5); NK-like (12); monocytes (10, 14); unclear (2, 6, 7, 11; cluster 11 likely doublets). No B cells (MS4A1 and CD79A absent).
- CARPos CD8 T cells are mostly in cluster 8; CARNeg mostly in clusters 3, 4, and 13. CAR status is confounded with cluster.
- Within-cluster CARNeg minus CARPos signature score: cluster 3 -0.20 (p = 0.005); cluster 8 +0.10 (p = 0.13); clusters 4 and 13 untestable (18 and 2 CARPos cells).
- Pooled comparison significant (p = 2.3e-5) but confounded by cluster. Not interpreted.
- Limits: one patient, so the statistical unit is the patient and p-values are descriptive. Gene set approximated from the abstract, not the authors' full signature.

## Open items

- [x] Confirm library type: 10x Genomics Next GEM 5′ v2 (GSM8248030)
- [x] Confirm Bag and Product are the same group: CAR T infusion product (GSM8248030)
- [x] Write the README and interview notes (step 5)
- [x] Cluster, annotate, and compare CARNeg with CARPos (exploratory, step 6)

## Environment

- **Python:** 3.12.4 in `.venv/`
- **Packages:** scanpy, anndata, pandas, duckdb
- **Machine:** MacBook Air, 8 GB RAM, 20 GB free storage
