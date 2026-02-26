# Report dashboard contract (Phase 6 Wave 2)

## Overview

Phase 6 Wave 2 adds a lightweight interactive dashboard package on top of report bundle artifacts.

Run with:

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host_variants.tsv> -profile standard`

## Inputs

Dashboard generation consumes Wave 1 bundle outputs:
- `report_bundle_metrics.json`
- `report_bundle_top_signals.tsv`

## Outputs

- `report_dashboard.html`
- `report_dashboard_data.json`

### `report_dashboard.html`

Self-contained interactive page with:
- metric cards for key WGS/RNA summary values
- sortable top-signal table
- source filter (`wgs_assoc`, `rna_assoc`, all)
- feature/gene text search

### `report_dashboard_data.json`

Machine-readable dashboard metadata:
- `metrics` (copied from report bundle)
- `signal_count`

## Notes

- Wave 2 dashboard is static/self-contained and requires no server runtime.
- Wave 3 addresses scheduler portability/deployment hardening.
