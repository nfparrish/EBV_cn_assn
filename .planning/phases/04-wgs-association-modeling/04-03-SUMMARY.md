---
phase: 04-wgs-association-modeling
plan: "03"
subsystem: association-diagnostics-qc
tags:
  - nextflow
  - association
  - diagnostics
  - qc
provides:
  - Reproducible diagnostics artifacts for EBV and host association tracks
  - Phase-level QC summary and top-hit table for sign-off
affects:
  - 04-01
  - 04-02
tech-stack:
  added:
    - Nextflow DSL2 modules: ASSOC_DIAGNOSTICS, ASSOC_QC_SUMMARY
    - Python scripts for diagnostics aggregation and QC synthesis
  patterns:
    - Pathway-level diagnostics rolled up into phase-level QC flags
key-files:
  created:
    - pipeline/modules/local/assoc_diagnostics/main.nf
    - pipeline/modules/local/assoc_qc_summary/main.nf
    - pipeline/scripts/build_assoc_diagnostics.py
    - pipeline/scripts/build_assoc_qc_summary.py
    - pipeline/docs/association-diagnostics.md
  modified:
    - pipeline/workflows/wgs_association.nf
key-decisions:
  - decision: Use deterministic coverage/significance proxies for diagnostics in v1 wave 3.
    rationale: Provides reproducible quality indicators before introducing richer graphical diagnostics.
  - decision: Emit high-level QC flags (`low_model_coverage_any_pathway`, `no_nominal_hits`).
    rationale: Creates actionable phase-signoff signals for small pilot and scale-up runs.
patterns-established:
  - Association workflow now emits execution stats, diagnostics, and QC summaries in one chain.
duration: 35min
completed: 2026-02-26
---

# Phase 4: WGS Association Modeling Summary (04-03)

Implemented Wave 3 diagnostics and QC reporting for association outputs.

## Accomplishments
- Added `ASSOC_DIAGNOSTICS` module to compute pathway-level modeled/unmodeled counts and summary metrics.
- Added `ASSOC_QC_SUMMARY` module to generate QC summary JSON, QC report markdown, and top-hits table.
- Wired diagnostics + QC modules into `WGS_ASSOCIATION` workflow emits.
- Added diagnostics documentation in `pipeline/docs/association-diagnostics.md`.

## Output contracts added
- `association_diagnostics.tsv`
- `association_diagnostics_summary.json`
- `association_qc_summary.json`
- `association_qc_report.md`
- `association_qc_top_hits.tsv`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow assoc --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants.tsv --host_variant_format tsv --outdir pipeline/tests/pilot_assoc_run_04_03 -profile standard`
- Diagnostics/QC modules completed successfully and emitted expected artifacts.
- Pilot QC correctly flags low coverage due to `association_min_samples=3` with only 2 manifest samples.

## Next Readiness
- Phase 4 is complete and ready for transition planning into Phase 5 RNA-seq integration.
