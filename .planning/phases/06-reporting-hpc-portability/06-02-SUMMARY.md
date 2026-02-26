---
phase: 06-reporting-hpc-portability
plan: "02"
subsystem: report-dashboard
tags:
  - reporting
  - dashboard
  - nextflow
  - html
provides:
  - Interactive dashboard package generated from report bundle outputs
  - Self-contained HTML with filter/sort/search over top signals
affects:
  - 06-01
tech-stack:
  added:
    - Nextflow DSL2 module: REPORT_DASHBOARD
    - Python renderer: build_report_dashboard.py
  patterns:
    - Workflow-native static dashboard packaging (no external server dependency)
key-files:
  created:
    - pipeline/modules/local/report_dashboard/main.nf
    - pipeline/scripts/build_report_dashboard.py
    - pipeline/docs/report-dashboard.md
    - .planning/phases/06-reporting-hpc-portability/06-02-PLAN.md
  modified:
    - pipeline/workflows/reporting_bundle.nf
key-decisions:
  - decision: Build a self-contained HTML dashboard artifact rather than a framework app in Wave 2.
    rationale: Keeps output portable and scheduler-friendly for batch-generated artifacts.
  - decision: Drive dashboard from Wave 1 report bundle contracts only.
    rationale: Preserves deterministic reporting boundaries and minimizes coupling.
patterns-established:
  - `--workflow report` now emits narrative bundle + interactive dashboard artifacts together.
duration: 35min
completed: 2026-02-26
---

# Phase 6: Reporting & HPC Portability Summary (06-02)

Implemented Wave 2 interactive dashboard packaging.

## Accomplishments
- Added `REPORT_DASHBOARD` module and integrated it into `REPORTING_BUNDLE`.
- Added dashboard renderer producing:
  - metric cards
  - sortable top-signals table
  - source filter
  - text search over feature/gene identifiers
- Added dashboard contract documentation.

## Output contracts added
- `report_dashboard.html`
- `report_dashboard_data.json`

## Validation Notes
- Successful end-to-end run:
  `nextflow run pipeline/main.nf --workflow report --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants_local.tsv --rna_input_mode matrix --rna_assoc_min_samples 2 --outdir pipeline/tests/pilot_report_run_06_02 -profile standard`
- Verified dashboard artifacts emitted and browser-renderable over local HTTP preview.

## Next Readiness
- Wave 3 (`06-03`) can finalize scheduler profiles, deployment docs, and reproducibility checklists.
