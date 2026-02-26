---
phase: 02-ebv-copy-number-engine
plan: "03"
subsystem: ebv-cn-confidence-qc
tags:
  - nextflow
  - ebv
  - cn
  - confidence
  - qc
provides:
  - Rule-based confidence score and tier per CN estimate
  - Cohort-level CN QC summary and low-confidence review table
affects:
  - 02-01
  - 02-02
tech-stack:
  added:
    - Nextflow DSL2 modules: EBV_CN_CONFIDENCE, EBV_CN_QC
    - Python helper scripts for confidence scoring and cohort QC reporting
  patterns:
    - Deterministic score-penalty confidence model with tiering
key-files:
  created:
    - pipeline/modules/local/ebv_cn_confidence/main.nf
    - pipeline/modules/local/ebv_cn_qc/main.nf
    - pipeline/scripts/build_ebv_cn_confidence.py
    - pipeline/scripts/build_ebv_cn_qc.py
    - pipeline/docs/ebv-cn-confidence.md
  modified:
    - pipeline/workflows/ebv_cn_engine.nf
key-decisions:
  - decision: Keep Wave 2 confidence model rule-based with explicit penalty thresholds.
    rationale: Matches user-selected method and keeps pilot behavior transparent and reproducible.
  - decision: Separate per-sample confidence annotation from cohort QC reporting.
    rationale: Preserves modularity and enables later replacement of confidence engine without changing QC outputs.
patterns-established:
  - Normalized CN outputs now feed confidence annotation, which then feeds pilot QC summary generation.
duration: 35min
completed: 2026-02-26
---

# Phase 2: EBV Copy Number Engine Summary (02-03)

Implemented confidence annotation and CN QC summaries for Wave 2.

## Accomplishments
- Added `EBV_CN_CONFIDENCE` module producing `ebv_cn_confidence.tsv` and `ebv_cn_confidence_summary.json`.
- Added `EBV_CN_QC` module producing `ebv_cn_qc_summary.json`, `ebv_cn_qc_report.md`, and `ebv_cn_low_confidence.tsv`.
- Wired confidence and QC modules into `pipeline/workflows/ebv_cn_engine.nf` and exposed outputs through workflow emits.
- Documented rule-based thresholds and interpretation in `pipeline/docs/ebv-cn-confidence.md`.

## Validation Notes
- Pilot run completed successfully with confidence and QC stages:
  `nextflow run pipeline/main.nf --workflow cn --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --outdir pipeline/tests/pilot_cn_run_wave2 -profile standard -resume`
- `EBV_CN_CONFIDENCE` and `EBV_CN_QC` both executed successfully after path-handling fix.

## Next Readiness
- Phase 2 now includes end-to-end CN generation, confidence scoring, and pilot QC reporting for downstream association integration.
