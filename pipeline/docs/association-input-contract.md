# Association input contract (Phase 4 Wave 1)

## Overview

Phase 4 Wave 1 prepares model-ready association inputs by harmonizing covariates and z-scoring EBV CN outcome values.

Run with:

`nextflow run pipeline/main.nf --workflow assoc --manifest <manifest.csv> -profile standard`

Optional external covariates input:

`--association_covariates_input <path/to/covariates.tsv>`

## Covariate harmonization output

`association_covariates.tsv` columns:

- `sample_id`
- `dataset`
- `population`
- `sex`
- `pc1`..`pc5`
- `sequence_lab`
- `batch`
- `lab_batch`
- `covariate_mode`

Auto fallback policy:

- If `sequence_lab` and `batch` are available, both are retained.
- If either is missing, `lab_batch` is synthesized from available values.

`association_covariates_report.json` records source mode, overlap counts, and schema fields.

## CN preparation output

`association_input.tsv` columns:

- `sample_id`
- `cn_value`
- `cn_zscore`
- covariates (`dataset`, `population`, `sex`, `pc1`..`pc5`, `sequence_lab`, `batch`, `lab_batch`)

`association_cn_prepare_report.json` records:

- selected CN source column (`association_cn_column`)
- overlap sample count
- CN mean/std used for z-score

## Relevant params

- `association_covariates_input` (optional)
- `association_cn_column` (default `depth_ratio_cn_norm`)
- `association_model_family` (planning default for downstream wave)
- `association_model_link` (planning default for downstream wave)
- `association_multiple_testing` (planning default for downstream wave)
