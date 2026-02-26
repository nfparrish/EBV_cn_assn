---
phase: 05-rna-seq-integration
plan: "03"
subsystem: rna-variant-association-qc
tags:
  - nextflow
  - rna-seq
  - association
  - variant
  - qc
provides:
  - EBV variant→expression association outputs for dual normalization tracks
  - Integrated RNA association QC summary/report/top-hit artifacts
affects:
  - 05-01
  - 05-02
tech-stack:
  added:
    - Nextflow DSL2 modules: RNA_ASSOC_VARIANT, RNA_ASSOC_QC
    - Python scripts for variant-expression modeling and integrated RNA QC aggregation
  patterns:
    - Dual-track RNA association with shared inferential schema and QC rollup
key-files:
  created:
    - pipeline/modules/local/rna_assoc_variant/main.nf
    - pipeline/modules/local/rna_assoc_qc/main.nf
    - pipeline/scripts/run_rna_assoc_variant.py
    - pipeline/scripts/build_rna_assoc_qc.py
    - pipeline/docs/rna-association-qc.md
  modified:
    - pipeline/workflows/rna_integration.nf
key-decisions:
  - decision: Generate EBV variant inputs directly inside `RNA_INTEGRATION` via `EBV_REFERENCE_PREP` + `EBV_VARIANT_CALL`.
    rationale: Keeps RNA workflow self-contained for Wave 3 without requiring external host-variant ingestion parameters.
  - decision: Run variant→expression associations for both normalization tracks with the same covariate policy used in CN→expression.
    rationale: Preserves comparability across RNA normalization branches.
  - decision: Gate inferential statistics on positive residual degrees-of-freedom.
    rationale: Prevents spurious tiny-cohort p-values (e.g., n=2) from appearing as false positive signal.
patterns-established:
  - `--workflow rna` now emits CN-association, variant-association, and integrated RNA-association QC outputs.
duration: 55min
completed: 2026-02-26
---

# Phase 5: RNA-seq Integration Summary (05-03)

Implemented Wave 3 EBV variant→expression association and RNA integration QC summaries.

## Accomplishments
- Added `RNA_ASSOC_VARIANT` module and dual invocation in `RNA_INTEGRATION`:
  - `log2_tpm_plus_1`
  - `quantile_normalized`
- Added internal EBV variant generation path in RNA workflow (`EBV_REFERENCE_PREP`, `EBV_VARIANT_CALL`) and collection for variant association.
- Added `RNA_ASSOC_QC` module to aggregate CN/variant association outputs into a unified QC summary/report/top-hits set.
- Added `pipeline/docs/rna-association-qc.md` documenting variant association and QC contracts.
- Added small-cohort inferential safeguard in variant association execution (emit `nan` stats when DOF is insufficient).

## Output contracts added
- Variant association (per normalization):
  - `rna_assoc_variant_results.tsv`
  - `rna_assoc_variant_report.json`
- Integrated RNA association QC:
  - `rna_assoc_qc_summary.json`
  - `rna_assoc_qc_report.md`
  - `rna_assoc_qc_top_hits.tsv`
- Workflow emits:
  - `rna_assoc_variant_log2_results`
  - `rna_assoc_variant_log2_report`
  - `rna_assoc_variant_quantile_results`
  - `rna_assoc_variant_quantile_report`
  - `rna_assoc_qc_summary`
  - `rna_assoc_qc_report`
  - `rna_assoc_qc_top_hits`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow rna --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --rna_assoc_min_samples 2 --outdir pipeline/tests/pilot_rna_run_05_03b -profile standard`
- Verified both variant association branches and integrated QC process execute successfully.
- Confirmed tiny-cohort guard behavior: variant and CN inferential fields are `nan` when model DOF is insufficient, and QC reflects low-coverage model health flags.

## Phase 5 Status
- Wave 1 (`05-01`) complete.
- Wave 2 (`05-02`) complete.
- Wave 3 (`05-03`) complete.
- Phase 5 baseline implementation is now complete and ready for phase-level wrap-up artifacts.
