---
phase: 06-reporting-hpc-portability
plan: "03"
subsystem: portability-reproducibility
tags:
  - hpc
  - reproducibility
  - reporting
  - nextflow
provides:
  - Automated reproducibility verification artifacts for report outputs
  - Scheduler/deployment portability guidance for local, slurm, and pbs
affects:
  - 06-01
  - 06-02
tech-stack:
  added:
    - Nextflow DSL2 module: REPORT_REPRO_CHECK
    - Python reproducibility script: build_report_reproducibility.py
  patterns:
    - Hash-based report artifact validation with manifest/checksum/report triad
key-files:
  created:
    - pipeline/modules/local/report_repro_check/main.nf
    - pipeline/scripts/build_report_reproducibility.py
    - pipeline/docs/hpc-portability.md
    - pipeline/docs/reproducibility-checks.md
    - .planning/phases/06-reporting-hpc-portability/06-03-PLAN.md
  modified:
    - pipeline/workflows/reporting_bundle.nf
key-decisions:
  - decision: Integrate reproducibility checks into `--workflow report` instead of post-hoc manual scripts.
    rationale: Ensures every report run emits verifiable checksums by default.
  - decision: Keep portability guidance profile-centric using existing `standard`, `slurm`, and `pbs` definitions.
    rationale: Aligns with current config model while giving actionable deployment steps.
patterns-established:
  - `--workflow report` now emits bundle, dashboard, and reproducibility verification artifacts in one run.
duration: 30min
completed: 2026-02-26
---

# Phase 6: Reporting & HPC Portability Summary (06-03)

Implemented Wave 3 scheduler portability and reproducibility checks.

## Accomplishments
- Added `REPORT_REPRO_CHECK` module and integrated it into `REPORTING_BUNDLE`.
- Added reproducibility script generating:
  - checksum table (SHA-256)
  - machine-readable manifest
  - human-readable reproducibility report
- Added HPC/deployment guidance and reproducibility documentation.

## Output contracts added
- `report_reproducibility_manifest.json`
- `report_reproducibility_checksums.tsv`
- `report_reproducibility_report.md`

## Validation Notes
- Successful run:
  `nextflow run pipeline/main.nf --workflow report --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants_local.tsv --rna_input_mode matrix --rna_assoc_min_samples 2 --outdir pipeline/tests/pilot_report_run_06_03 -profile standard`
- Verified all reproducibility outputs are emitted with expected schemas and non-empty checksums.

## Phase 6 Status
- Wave 1 (`06-01`) complete.
- Wave 2 (`06-02`) complete.
- Wave 3 (`06-03`) complete.
- Reporting + HPC portability baseline is complete.
