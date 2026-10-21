# Interview notes: CAR-T single-cell QC project

## 30-second summary

I built a QC and metadata-validation pipeline on public single-cell data from a 2024 *Blood* CAR-T study. I analyzed one patient, 12,244 cells, and kept 72% after filtering on gene counts and mitochondrial reads. I validated the metadata against 15 project-defined checks, and I recorded every decision in a lineage log. I then ran an exploratory test of the paper's bystander T-cell signature. It didn't reproduce cleanly, and the reason was a confound between CAR status and cell state.

## Decisions and why

- **One patient (Pt011):** It was the only patient with CAR-positive, CAR-negative, and PBMC groups together. Pt006 looked complete at first, but I found it had no CAR-negative group and switched.
- **15% mitochondrial cutoff:** I tested 10%, 15%, and 20%. At 10%, I'd lose 20% of CAR-positive cells. At 20%, I'd keep more likely-damaged cells. 15% kept 88.7% of CAR-positive cells.
- **The cost:** The filter removed 43% of the infusion product cells. Their higher mitochondrial reads may be real biology, so I flagged it rather than calling it noise.
- **One alias:** The source labels the infusion product "Bag" in one column and "Product" in another. I confirmed in the GEO sample record that they're the same group, and I documented it in the validation report.

## The signature check

- Within matched clusters, results were mixed. Cluster 3 went the opposite way from the paper (p = 0.005). Cluster 8 showed a small, non-significant increase (p = 0.13).
- The pooled result was significant (p = 2×10⁻⁵), but I didn't report it. CAR-positive cells sit mostly in one cluster, so the pooled test measured cell state, not CAR status.
- The gene set came from the paper's abstract, not its full signature, and the study has one patient. So the result is descriptive, not a test of the paper's claim.

## Likely questions

**Why did you filter this way?** Thresholds were chosen by looking at how each cell group was affected, not by overall retention alone.

**Did you reproduce the paper?** No. The signature check was exploratory, on one patient, with an approximate gene set. The within-cluster results didn't support the claim.

**What went wrong?** My first pass said Pt006 had all four groups. I checked the counts and found it didn't. I switched to Pt011, and I corrected the earlier claim in the README and lineage record.

**Why did the pooled result look significant but mean nothing?** Confounding. The CAR-positive cells came mostly from one cell cluster, and the scores differ by cluster more than by CAR status. Comparing within clusters is the correct test, and that's where the signal disappeared.

**What's the statistical unit?** The patient, which is one here. The cell-level p-values assume the cells are independent, which they aren't, so they're descriptive only.

**What would you do next?** Add patients with CAR-negative groups, use the paper's full signature, match cell states before comparing, and look at the TCR data, which covers about half the cells.

**What didn't you do?** No raw reads (dbGaP submission pending), no differential expression, no formal cell-type annotation beyond marker genes. Several clusters were unclear, and cluster 11 probably contains doublets.

## What to say about the tools

I used scanpy, AnnData, pandas, and Leiden clustering. I did not use dbt, Snowflake, or Airflow on this project. Those are in the NYC Subway project.
