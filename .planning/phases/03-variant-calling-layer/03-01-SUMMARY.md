---
phase: 03-variant-calling-layer
plan: "01"
subsystem: ebv-variant-call-matrix
tags:
  - nextflow
  - ebv
  - variants
  - mutect2
  - matrix
provides:
  - Per-sample EBV variant call artifacts with provenance fields
  - Cohort AF and binary EBV variant matrices
affects: []
tech-stack:
  added:
    - Nextflow DSL2 workflow: EBV_VARIANT_LAYER
    - Nextflow DSL2 modules: EBV_VARIANT_CALL, EBV_VARIANT_MATRIX
    - Python scripts for per-sample variant emission and cohort matrix synthesis
  patterns:
    - Call-per-sample then collect-and-aggregate matrix strategy
key-files:
  created:
    - pipeline/workflows/ebv_variant_layer.nf
    - pipeline/modules/local/ebv_variant_call/main.nf
    - pipeline/modules/local/ebv_variant_matrix/main.nf
    - pipeline/scripts/compute_ebv_variants_mutect2.py
    - pipeline/scripts/build_ebv_variant_matrices.py
    - pipeline/docs/ebv-variant-calling.md
  modified:
    - pipeline/main.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Implement Wave 1 EBV calling around Mutect2-baseline provenance tags.
    rationale: Aligns with user-selected Phase 3 caller strategy while keeping output contract deterministic for pilot progression.
  - decision: Emit both AF and binary matrices in the same aggregation pass.
    rationale: Satisfies dual-output requirement and avoids divergence between derived artifacts.
patterns-established:
  - `--workflow variant` route extends main entrypoint without impacting intake/CN paths.
duration: 40min
completed: 2026-02-26
---

# Phase 3: Variant Calling Layer Summary (03-01)

Implemented Wave 1 EBV variant calling and cohort matrix generation.

## Accomplishments
- Added new workflow path `--workflow variant` in `pipeline/main.nf`.
- Added `EBV_VARIANT_LAYER` workflow to orchestrate reference prep, per-sample variant calling, and cohort matrix aggregation.
- Added `EBV_VARIANT_CALL` module emitting per-sample `*.ebv_variants.tsv` artifacts with variant and metadata/provenance fields.
- Added `EBV_VARIANT_MATRIX` module emitting:
  - `ebv_variant_catalog.tsv`
  - `ebv_variant_af_matrix.tsv`
  - `ebv_variant_binary_matrix.tsv`
  - `ebv_variant_matrix_summary.json`
- Added config knobs in `pipeline/conf/base.config` for variant caller and binary threshold defaults.
- Documented outputs and behavior in `pipeline/docs/ebv-variant-calling.md`.

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow variant --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --outdir pipeline/tests/pilot_variant_run -profile standard`
- `EBV_VARIANT_CALL` and `EBV_VARIANT_MATRIX` completed successfully for all pilot samples.

## Next Readiness
- Phase 3 Wave 2 (`03-02`) can now integrate host variant ingestion/harmonization into this variant layer.
