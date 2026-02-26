# Association diagnostics and QC (Phase 4 Wave 3)

## Overview

Phase 4 Wave 3 adds model diagnostics and QC summarization on top of Wave 2 association outputs.

Added modules:

- `ASSOC_DIAGNOSTICS`
- `ASSOC_QC_SUMMARY`

These consume:

- `assoc_ebv_results.tsv`
- `assoc_host_results.tsv`

## Diagnostics outputs

- `association_diagnostics.tsv`
- `association_diagnostics_summary.json`

Diagnostics include per-pathway:

- total/modelled/unmodelled variant counts
- median absolute effect size and median standard error
- median absolute z-score
- nominal (`p<0.05`) and FDR (`<0.10`) counts
- model coverage proxy rates

## QC outputs

- `association_qc_summary.json`
- `association_qc_report.md`
- `association_qc_top_hits.tsv`

QC summary highlights:

- pathway-level model coverage and signal counts
- top hits table for quick inspection
- high-level flags (e.g., low coverage / no nominal hits)

## Notes

Wave 3 diagnostics are intentionally lightweight and reproducible for phase sign-off.
They are designed to be extended with richer residual plots and calibration checks in future iterations.
