---
phase: 04-wgs-association-modeling
plan: "01"
subsystem: association-input-prep
tags:
  - nextflow
  - association
  - covariates
  - cn
  - zscore
provides:
  - Harmonized association covariate matrix with lab/batch auto fallback
  - Z-scored CN association input table and transform metadata
affects: []
tech-stack:
  added:
    - Nextflow DSL2 workflow: WGS_ASSOCIATION
    - Nextflow DSL2 modules: ASSOC_COVARIATE_HARMONIZE, ASSOC_CN_PREPARE
    - Python scripts for covariate harmonization and CN z-score prep
  patterns:
    - Reuse Phase 2 CN outputs as upstream inputs to model-ready preparation
key-files:
  created:
    - pipeline/workflows/wgs_association.nf
    - pipeline/modules/local/assoc_covariate_harmonize/main.nf
    - pipeline/modules/local/assoc_cn_prepare/main.nf
    - pipeline/scripts/build_assoc_covariates.py
    - pipeline/scripts/prepare_assoc_cn.py
    - pipeline/docs/association-input-contract.md
  modified:
    - pipeline/main.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Introduce `--workflow assoc` as the Phase 4 entrypoint while reusing `EBV_CN_ENGINE` outputs.
    rationale: Keeps Wave 1 focused on model-input readiness without duplicating CN logic.
  - decision: Default CN source for z-score prep is `depth_ratio_cn_norm` via `association_cn_column`.
    rationale: Provides deterministic baseline while allowing configurable source column overrides.
patterns-established:
  - Association prep now emits reproducible covariate and z-score artifacts ready for GLM fitting in Wave 2.
duration: 45min
completed: 2026-02-26
---

# Phase 4: WGS Association Modeling Summary (04-01)

Implemented Wave 1 covariate harmonization and CN z-score preparation.

## Accomplishments
- Added new `assoc` workflow route and `WGS_ASSOCIATION` orchestration.
- Added `ASSOC_COVARIATE_HARMONIZE` module and `association_covariates.tsv` + report outputs.
- Added `ASSOC_CN_PREPARE` module and `association_input.tsv` + z-score metadata report outputs.
- Added Phase 4 config params for covariate input and association model defaults.
- Documented contracts in `pipeline/docs/association-input-contract.md`.

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow assoc --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --outdir pipeline/tests/pilot_assoc_run -profile standard -resume`
- `ASSOC_COVARIATE_HARMONIZE` and `ASSOC_CN_PREPARE` completed successfully after staged-path fix.

## Next Readiness
- Wave 2 (`04-02`) can now implement EBV and host variant GLM execution using this model-ready input layer.
