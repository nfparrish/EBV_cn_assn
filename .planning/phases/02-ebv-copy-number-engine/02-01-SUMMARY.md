---
phase: 02-ebv-copy-number-engine
plan: "01"
subsystem: ebv-reference-depth
tags:
  - nextflow
  - ebv
  - cn
  - reference
  - depth
provides:
  - EBV reference manifest contract
  - Per-sample EBV depth intermediate artifacts
affects: []
tech-stack:
  added:
    - Nextflow DSL2 modules: EBV_REFERENCE_PREP, EBV_DEPTH_EXTRACT
    - Python helper scripts for deterministic reference/depth artifacts
  patterns:
    - Reference-manifest-driven sample metric extraction
key-files:
  created:
    - pipeline/modules/local/ebv_reference_prep/main.nf
    - pipeline/modules/local/ebv_depth_extract/main.nf
    - pipeline/scripts/build_ebv_reference_manifest.py
    - pipeline/scripts/compute_ebv_depth.py
    - pipeline/docs/ebv-reference.md
  modified:
    - pipeline/workflows/ebv_cn_engine.nf
    - pipeline/main.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Support both default accession and optional user FASTA modes through one reference manifest artifact.
    rationale: Keeps reference provenance explicit and portable across compute environments.
  - decision: Emit one deterministic depth JSON per sample as the canonical Wave 1 intermediate.
    rationale: Simplifies downstream CN computation and troubleshooting.
patterns-established:
  - CN workflow is selectable via `--workflow cn` without affecting intake workflow behavior.
duration: 40min
completed: 2026-02-26
---

# Phase 2: EBV Copy Number Engine Summary (02-01)

Implemented EBV reference prep and per-sample depth extraction foundations for the CN workflow path.

## Accomplishments
- Added `EBV_REFERENCE_PREP` to emit `ebv_reference_manifest.json` with explicit source metadata.
- Added `EBV_DEPTH_EXTRACT` to generate per-sample `*.ebv_depth.json` artifacts with required metadata and support metrics.
- Added workflow routing in `pipeline/main.nf` so CN execution is selected via `--workflow cn`.
- Added reference configuration knobs in `pipeline/conf/base.config` (`reference_accession`, `viral_reference_fasta`, baseline/method defaults).
- Documented reference behavior in `pipeline/docs/ebv-reference.md`.

## Validation Notes
- Ran local pilot: `nextflow run pipeline/main.nf --workflow cn --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --outdir pipeline/tests/pilot_cn_run -profile standard`.
- `EBV_REFERENCE_PREP` and `EBV_DEPTH_EXTRACT` completed successfully for all pilot samples.

## Next Readiness
- Inputs required for CN estimation are now present and wired for Wave 1 module expansion (`02-02`).
