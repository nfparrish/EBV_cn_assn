# RNA variant presence QC contract

This module validates expected EBV coding variants and donor-specific coding variants per RNA donor/sample.

## Process
- `RNA_VARIANT_PRESENCE_QC`

## Inputs
- Per-sample EBV variant TSV files from `EBV_VARIANT_CALL`.
- `rna_sample_map.tsv` from RNA ingest.
- Optional: `--rna_expected_donor_variants_input` TSV with columns:
  - `sample_id`
  - `variant_id`

## Config parameters
- `rna_expected_coding_variants` (comma-separated variant IDs)
- `rna_expected_donor_variants_input` (optional TSV)
- Uses variant thresholds:
  - `variant_filter_ebv_min_af`
  - `variant_filter_ebv_min_dp`
  - `variant_filter_ebv_min_qual`

## Outputs

### `rna_variant_presence_checks.tsv`
Columns:
- `sample_id`
- `check_class` (`ebv_expected_coding` or `donor_specific_coding`)
- `variant_id`
- `present`
- `evidence_af`
- `evidence_dp`
- `evidence_qual`

### `rna_variant_presence_summary.json`
Contains:
- sample-level totals/present counts for both check classes
- filter thresholds used
- donor expectation source (`input_table` or `inferred_unique_variants`)
- global pass/fail flags

### `rna_variant_presence_report.md`
Human-readable report of scope and global flags.
