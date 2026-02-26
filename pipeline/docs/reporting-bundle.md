# Reporting bundle contract (Phase 6 Wave 1)

## Overview

Phase 6 Wave 1 adds a reporting workflow that aggregates WGS association and RNA integration outputs into publication-oriented bundle artifacts.

Run with:

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host_variants.tsv> -profile standard`

## Inputs

The reporting workflow executes and consumes:
- `WGS_ASSOCIATION` outputs (diagnostics, QC summary, top hits)
- `RNA_INTEGRATION` outputs (RNA QC summary, top hits, variant presence summary, ingest report)

## Outputs

- `report_bundle_summary.md`
- `report_bundle_metrics.json`
- `report_bundle_top_signals.tsv`

### `report_bundle_metrics.json`

Includes summary metrics grouped by:
- `wgs_assoc`
- `rna_assoc`
- `bundle`

### `report_bundle_top_signals.tsv`

Merged top signals across WGS and RNA layers with fields:
- `source`
- `pathway`
- `layer`
- `feature_id`
- `gene_id`
- `effect_size`
- `p_value`
- `p_value_fdr`
- `n_samples`

## Notes

- Wave 1 focuses on deterministic markdown/json/tsv bundles.
- Interactive dashboard packaging and scheduler/deployment hardening are handled in later Phase 6 waves.
