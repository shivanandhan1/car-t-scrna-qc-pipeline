# Metadata validation: Pt011 filtered data

| Check | Result | Detail |
|---|---|---|
| sample_id present and complete | PASS | 0 missing |
| donor present and complete | PASS | 0 missing |
| origin present and complete | PASS | 0 missing |
| Patient present and complete | PASS | 0 missing |
| Cell_Source present and complete | PASS | 0 missing |
| Formatted_Sample_Name present and complete | PASS | 0 missing |
| gem_well present and complete | PASS | 0 missing |
| Cell_Source uses expected labels | PASS | all expected |
| Single patient (Pt011) | PASS | found: ['Pt011'] |
| donor matches Patient | PASS | 0 mismatches |
| sample_id = Patient_Cell_Source | PASS | 0 mismatches |
| Formatted_Sample_Name matches sample_id (alias applied) | PASS | 0 mismatches after alias; alias used for 4228 cells |
| cell barcodes unique | PASS | 0 duplicates |
| cell barcode format (e.g. AAACCTGAGGAATTAC-1) | PASS | 0 malformed |
| gem_well is a positive integer | PASS | 0 invalid |
