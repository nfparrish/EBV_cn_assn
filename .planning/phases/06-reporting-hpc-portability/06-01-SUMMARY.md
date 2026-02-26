---
phase: 06-reporting-hpc-portability
plan: "01"
subsystem: reporting-bundle
tags:
  - reporting
  - nextflow
  - markdown
  - qc
  - publication
provides:
  - Workflow-native reporting bundle artifacts from WGS association and RNA integration outputs
  - Merged top-signal table and machine-readable summary metrics
affects:
  - 04-03
  - 05-03
tech-stack:
  added:
    - Nextflow DSL2 workflow: REPORTING_BUNDLE
    - Nextflow DSL2 module: REPORT_BUNDLE
    - Python report builder: build_report_bundle.py
  patterns:
    - Deterministic markdown/json/tsv reporting synthesis from existing workflow contracts
key-files:
  created:
    - pipeline/workflows/reporting_bundle.nf
    - pipeline/modules/local/report_bundle/main.nf
    - pipeline/scripts/build_report_bundle.py
    - pipeline/docs/reporting-bundle.md
    - .planning/phases/06-reporting-hpc-portability/06-CONTEXT.md
    - .planning/phases/06-reporting-hpc-portability/06-RESEARCH.md
    - .planning/phases/06-reporting-hpc-portability/06-01-PLAN.md
  modified:
    - pipeline/main.nf
key-decisions:
  - decision: Add `--workflow report` as a first-class entrypoint rather than external post-processing.
    rationale: Keeps reporting reproducible and scheduler-portable under the same workflow framework.
  - decision: Build bundle from stable upstream contracts (`assoc` and `rna` outputs) without introducing new statistical methods.
    rationale: Minimizes methodological drift and keeps Wave 1 focused on packaging/synthesis.
patterns-established:
  - Reporting bundle execution is now workflow-native and can be generated directly from a manifest-driven run.
duration: 40min
completed: 2026-02-26
---

# Phase 6: Reporting & HPC Portability Summary (06-01)

Implemented Wave 1 reporting bundle generation.

## Accomplishments
- Added `REPORTING_BUNDLE` workflow and `--workflow report` router support.
- Added `REPORT_BUNDLE` module to synthesize WGS/RNA QC artifacts into:
  - narrative markdown summary
  - machine-readable metrics JSON
  - merged top-signals TSV
- Added reporting contract documentation in `pipeline/docs/reporting-bundle.md`.
- Added Phase 6 planning context/research/plan files for continuity.

## Output contracts added
- `report_bundle_summary.md`
- `report_bundle_metrics.json`
- `report_bundle_top_signals.tsv`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow report --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants_local.tsv --rna_input_mode matrix --rna_assoc_min_samples 2 --outdir pipeline/tests/pilot_report_run_06_01 -profile standard`
- Verified all report bundle artifacts are emitted with expected schema.

## Next Readiness
- Wave 2 (`06-02`) can implement interactive dashboard packaging from report bundle outputs.
