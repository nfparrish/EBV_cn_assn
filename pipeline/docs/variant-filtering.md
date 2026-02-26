# Variant filtering and provenance (Phase 3 Wave 3)

## Overview

Phase 3 Wave 3 adds configurable balanced filtering and explicit provenance outputs across EBV and host variant pathways.

Workflow additions:

- `VARIANT_FILTER`
- `VARIANT_PROVENANCE`

These modules run inside `--workflow variant` after EBV matrix generation and host harmonization.

## Configurable thresholds

Defaults (balanced):

- `variant_filter_profile = balanced`
- `variant_filter_ebv_min_af = 0.05`
- `variant_filter_ebv_min_dp = 20`
- `variant_filter_ebv_min_qual = 60`
- `variant_filter_host_min_dosage = 1.0`

All thresholds can be overridden at runtime via `--param value`.

## Filter outputs

- `ebv_variants_unfiltered.tsv`
- `ebv_variants_filtered.tsv`
- `host_variant_filtered.tsv`
- `variant_filter_rules.json`
- `variant_filter_counts.json`

## Provenance outputs

- `variant_provenance_summary.json`
- `variant_provenance_counts.tsv`
- `variant_provenance_report.md`

These outputs record retained/excluded counts and sample-overlap context for VAR-03 reproducibility.

## Interpretation

- `variant_filter_counts.json` is the machine-readable source of retained/excluded counts by pathway.
- `variant_provenance_counts.tsv` is a tabular ledger suitable for audit checks and quick diffs.
- `variant_provenance_report.md` is human-facing phase-signoff text for pilot validation.
