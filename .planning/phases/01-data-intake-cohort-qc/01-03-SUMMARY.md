---
phase: 01-data-intake-cohort-qc
plan: "03"
subsystem: qc-reporting
tags:
  - qc
  - multiqc
  - cohort-aggregation
provides:
  - Per-sample QC metric artifacts
  - Cohort-level QC summary outputs and report contract
affects:
  - 01-01
  - 01-02
tech-stack:
  added:
    - QC_SAMPLE module
    - QC_COHORT module
    - MultiQC configuration for cohort table ingestion
  patterns:
    - Warning-first QC policy (no hard exclusions in phase 1)
key-files:
  created:
    - pipeline/modules/local/qc_sample/main.nf
    - pipeline/modules/local/qc_cohort/main.nf
    - pipeline/assets/multiqc/multiqc_config.yaml
    - pipeline/docs/qc-reporting.md
  modified:
    - pipeline/workflows/wgs_intake.nf
key-decisions:
  - decision: Keep phase-1 QC in warning-first mode.
    rationale: Allows pilot calibration before introducing hard exclusion thresholds.
  - decision: Produce both machine-readable tables and narrative report markdown.
    rationale: Supports downstream automation and human review.
patterns-established:
  - QC artifacts are generated per sample and aggregated at cohort level in the same workflow pass.
duration: 50min
completed: 2026-02-26
---

# Phase 1: Data Intake & Cohort QC Summary

Implemented Wave 2 QC metric generation and cohort reporting pipeline outputs.

## Performance
- **Duration:** 50min
- **Tasks:** 2 completed
- **Files modified:** 5

## Accomplishments
- Added `QC_SAMPLE` module producing structured per-sample metrics (`*.qc.json`).
- Added `QC_COHORT` module producing `qc_cohort_metrics.tsv`, `qc_cohort_summary.json`, and `qc_cohort_report.md`.
- Wired QC generation/aggregation into `WGS_INTAKE` workflow emits.
- Added MultiQC configuration and explicit QC output contract documentation.

## Validation Notes
- Nextflow binary is not available in this environment, so direct workflow execution was not possible.
- A deterministic local pilot fixture run validated expected output artifact contract under `pipeline/tests/pilot_outputs/`.
- `qc_cohort_summary.json` syntax validation passed.

## Files Created/Modified
- `pipeline/modules/local/qc_sample/main.nf`
- `pipeline/modules/local/qc_cohort/main.nf`
- `pipeline/assets/multiqc/multiqc_config.yaml`
- `pipeline/docs/qc-reporting.md`
- `pipeline/workflows/wgs_intake.nf`

## Next Phase Readiness
- Phase 1 now has intake, resumability, and QC reporting foundations complete.
- Phase transition/verification can proceed before moving into EBV copy-number engine implementation.
