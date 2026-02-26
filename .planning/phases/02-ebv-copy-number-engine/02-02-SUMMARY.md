---
phase: 02-ebv-copy-number-engine
plan: "02"
subsystem: ebv-cn-compute-normalize
tags:
  - nextflow
  - ebv
  - cn
  - normalization
provides:
  - Per-sample dual-method CN estimates with formula trace
  - Cohort raw and normalized CN tables + summary
affects:
  - 02-01
tech-stack:
  added:
    - Nextflow DSL2 modules: EBV_CN_COMPUTE, EBV_CN_NORMALIZE
    - Python helper scripts for CN formula and cohort normalization
  patterns:
    - Dual-method CN intermediates normalized at cohort level
key-files:
  created:
    - pipeline/modules/local/ebv_cn_compute/main.nf
    - pipeline/modules/local/ebv_cn_normalize/main.nf
    - pipeline/scripts/compute_ebv_cn_raw.py
    - pipeline/scripts/normalize_ebv_cn.py
    - pipeline/docs/ebv-cn-output-contract.md
  modified:
    - pipeline/workflows/ebv_cn_engine.nf
key-decisions:
  - decision: Implement both selected CN methods (`depth_ratio_cn` and `read_count_ratio_cn`) in all raw and normalized outputs.
    rationale: Aligns with user-selected method strategy and enables cross-method confidence in Wave 2.
  - decision: Normalize each method by cohort mean in Wave 1.
    rationale: Provides stable comparable features for confidence/QC expansion before introducing covariate-aware normalization.
patterns-established:
  - CN raw JSON and normalized cohort TSV artifacts form the handoff contract to confidence scoring.
duration: 45min
completed: 2026-02-26
---

# Phase 2: EBV Copy Number Engine Summary (02-02)

Implemented dual-method CN computation and cohort normalization outputs.

## Accomplishments
- Added `EBV_CN_COMPUTE` to produce per-sample `*.ebv_cn_raw.json` with formula traces and support metrics.
- Added `EBV_CN_NORMALIZE` to emit `ebv_cn_raw_cohort.tsv`, `ebv_cn_normalized.tsv`, and `ebv_cn_cohort_summary.json`.
- Wired CN module chain in `pipeline/workflows/ebv_cn_engine.nf` from reference prep → depth extraction → CN compute → cohort normalization.
- Documented output schema and contracts in `pipeline/docs/ebv-cn-output-contract.md`.

## Validation Notes
- CN pilot workflow completed successfully end-to-end with all Wave 1 modules.
- Raw and normalized CN artifacts were produced in work directories, confirming module wiring and file contracts.

## Next Readiness
- Wave 2 (`02-03`) can now implement rule-based confidence tiers and CN QC summaries on top of normalized CN outputs.
