---
phase: 05-rna-seq-integration
plan: "01"
subsystem: rna-ingest-normalize
tags:
  - nextflow
  - rna-seq
  - normalization
  - matrix
provides:
  - Canonical RNA expression + sample-map ingestion outputs
  - Dual normalization outputs: log2(TPM+1) and quantile-normalized matrices
affects: []
tech-stack:
  added:
    - Nextflow DSL2 workflow: RNA_INTEGRATION
    - Nextflow DSL2 modules: RNA_INGEST, RNA_NORMALIZE
    - Python scripts for ingestion and dual normalization
  patterns:
    - Canonical matrix first, parallel normalization branch outputs
key-files:
  created:
    - pipeline/workflows/rna_integration.nf
    - pipeline/modules/local/rna_ingest/main.nf
    - pipeline/modules/local/rna_normalize/main.nf
    - pipeline/scripts/build_rna_ingest.py
    - pipeline/scripts/normalize_rna_expression.py
    - pipeline/docs/rna-input-normalization.md
  modified:
    - pipeline/main.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Add `--workflow rna` as dedicated Phase 5 entrypoint.
    rationale: Keeps RNA integration isolated from WGS-only paths and easier to validate iteratively.
  - decision: Support optional external RNA matrix/metadata with deterministic fallback generation.
    rationale: Enables immediate pilot validation while preserving user-supplied data path for production runs.
patterns-established:
  - RNA ingestion now emits canonical matrix + sample map and two normalization artifacts with metadata.
duration: 40min
completed: 2026-02-26
---

# Phase 5: RNA-seq Integration Summary (05-01)

Implemented Wave 1 RNA-seq ingestion and dual normalization baseline.

## Accomplishments
- Added `RNA_INTEGRATION` workflow and `--workflow rna` route in `pipeline/main.nf`.
- Added `RNA_INGEST` module producing canonical expression matrix, sample map, and ingest report.
- Added `RNA_NORMALIZE` module producing both `log2(TPM+1)` and quantile-normalized matrices.
- Added Phase 5 RNA params in config and documented contracts in `pipeline/docs/rna-input-normalization.md`.

## Output contracts added
- `rna_expression_canonical.tsv`
- `rna_sample_map.tsv`
- `rna_ingest_report.json`
- `rna_expression_log2tpm.tsv`
- `rna_expression_quantile.tsv`
- `rna_normalization_report.json`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow rna --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --outdir pipeline/tests/pilot_rna_run_05_01 -profile standard`
- `RNA_INGEST` and `RNA_NORMALIZE` both completed successfully.

## Next Readiness
- Wave 2 (`05-02`) can now implement CN→expression association models over both normalization outputs.
