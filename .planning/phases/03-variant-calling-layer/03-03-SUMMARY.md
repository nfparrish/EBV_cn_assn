---
phase: 03-variant-calling-layer
plan: "03"
subsystem: variant-filtering-provenance
tags:
  - nextflow
  - variants
  - filtering
  - provenance
  - audit
provides:
  - Configurable balanced filtering across EBV and host pathways
  - Machine-readable and human-readable retained/excluded provenance artifacts
affects:
  - 03-01
  - 03-02
tech-stack:
  added:
    - Nextflow DSL2 modules: VARIANT_FILTER, VARIANT_PROVENANCE
    - Python scripts for filtering and provenance synthesis
  patterns:
    - Explicit threshold snapshot + stage-level count ledger
key-files:
  created:
    - pipeline/modules/local/variant_filter/main.nf
    - pipeline/modules/local/variant_provenance/main.nf
    - pipeline/scripts/filter_variants.py
    - pipeline/scripts/build_variant_provenance.py
    - pipeline/docs/variant-filtering.md
  modified:
    - pipeline/workflows/ebv_variant_layer.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Apply balanced defaults while exposing all threshold values as runtime params.
    rationale: Satisfies reproducibility and calibration needs without hard-coding policy.
  - decision: Keep provenance as first-class outputs (`json`, `tsv`, `md`).
    rationale: Supports both automated audits and human phase-signoff review.
patterns-established:
  - Variant workflow now emits pre/post-filter datasets and per-stage lineage counts for both EBV and host paths.
duration: 35min
completed: 2026-02-26
---

# Phase 3: Variant Calling Layer Summary (03-03)

Implemented Wave 3 balanced filtering and provenance reporting.

## Accomplishments
- Added `VARIANT_FILTER` module to apply configurable EBV and host thresholds.
- Added `VARIANT_PROVENANCE` module to emit stage/rule retained-excluded counts and narrative report.
- Wired filtering/provenance modules into `EBV_VARIANT_LAYER` and expanded workflow emits.
- Added filter threshold config parameters and operator docs.

## Output contracts added
- `ebv_variants_unfiltered.tsv`
- `ebv_variants_filtered.tsv`
- `host_variant_filtered.tsv`
- `variant_filter_rules.json`
- `variant_filter_counts.json`
- `variant_provenance_summary.json`
- `variant_provenance_counts.tsv`
- `variant_provenance_report.md`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow variant --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants.tsv --host_variant_format tsv --outdir pipeline/tests/pilot_variant_run_03_03 -profile standard`
- New `VARIANT_FILTER` and `VARIANT_PROVENANCE` modules completed successfully with EBV + host paths.

## Next Readiness
- Phase 3 now satisfies EBV matrix generation, host harmonization, and reproducible filtering/provenance requirements for transition into Phase 4 association modeling.
