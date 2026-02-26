---
phase: 05-rna-seq-integration
plan: "02"
subsystem: rna-cn-association
tags:
  - nextflow
  - rna-seq
  - association
  - cn
  - fdr
provides:
  - CN→expression association outputs for both normalization tracks
  - Standardized inferential schema and per-track reports
affects:
  - 05-01
tech-stack:
  added:
    - Nextflow DSL2 module: RNA_ASSOC_CN
    - Python script for CN-expression linear association execution with BH-FDR
  patterns:
    - Dual-track model invocation with shared output schema
key-files:
  created:
    - pipeline/modules/local/rna_assoc_cn/main.nf
    - pipeline/scripts/run_rna_assoc_cn.py
    - pipeline/docs/rna-association-contract.md
  modified:
    - pipeline/workflows/rna_integration.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Execute CN-expression associations for both log2(TPM+1) and quantile-normalized matrices in one workflow run.
    rationale: Ensures comparable outputs across normalization tracks and avoids divergent execution paths.
  - decision: Add adaptive covariate complexity for very small pilot cohorts.
    rationale: Prevents over-parameterized covariate matrices from blocking pilot execution while preserving richer covariates at scale.
patterns-established:
  - `--workflow rna` now includes ingestion, normalization, and CN-expression association outputs with per-track reports.
duration: 45min
completed: 2026-02-26
---

# Phase 5: RNA-seq Integration Summary (05-02)

Implemented Wave 2 CN→expression association modeling across both normalization tracks.

## Accomplishments
- Added `RNA_ASSOC_CN` module and dual invocation in `RNA_INTEGRATION` for:
  - `log2_tpm_plus_1`
  - `quantile_normalized`
- Added CN-expression association script with standardized result schema and BH-FDR adjustment fields.
- Added RNA association parameters in config (`rna_assoc_model_family`, `rna_assoc_model_link`, `rna_assoc_multiple_testing`, `rna_assoc_min_samples`).
- Added output contract documentation in `pipeline/docs/rna-association-contract.md`.

## Output contracts added
- `rna_assoc_cn_results.tsv` (per normalization track)
- `rna_assoc_cn_report.json` (per normalization track)
- Workflow emits:
  - `rna_assoc_cn_log2_results`
  - `rna_assoc_cn_log2_report`
  - `rna_assoc_cn_quantile_results`
  - `rna_assoc_cn_quantile_report`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow rna --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --rna_assoc_min_samples 2 --outdir pipeline/tests/pilot_rna_run_05_02b -profile standard`
- Both CN-association branches completed successfully.
- Pilot reports correctly indicate `min_p_value: null` when finite association statistics are not available in very small overlap sets.

## Next Readiness
- Wave 3 (`05-03`) can now implement EBV variant→expression association and RNA integration QC summaries.
