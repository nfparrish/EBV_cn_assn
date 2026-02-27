# RNA association output contract (Phase 5 Wave 2)

## Overview

Phase 5 Wave 2 adds CN→expression association modeling across both RNA normalization tracks.

Workflow module:

- `RNA_ASSOC_CN`

Invoked for:

- `log2_tpm_plus_1`
- `quantile_normalized`

## Required inputs

- `cn_normalized_table` from Phase 2 CN workflow
- `rna_expression_log2tpm.tsv` and `rna_expression_quantile.tsv`
- `rna_sample_map.tsv`

## Result outputs

- `rna_assoc_cn_results.tsv` (per invocation)
- `rna_assoc_cn_report.json` (per invocation)

Within `RNA_INTEGRATION` emits:

- `rna_assoc_cn_log2_results`
- `rna_assoc_cn_log2_report`
- `rna_assoc_cn_quantile_results`
- `rna_assoc_cn_quantile_report`

## Result schema

Columns in result TSV:

- `normalization`
- `gene_id`
- `n_samples`
- `effect_size`
- `std_error`
- `t_stat`
- `p_value`
- `p_value_fdr`
- `model_family`
- `model_link`
- `covariate_terms`

## Model baseline

- Linear model baseline with gaussian/identity defaults.
- Multiple testing uses Benjamini-Hochberg FDR.
- Covariate set includes empirically available PCs, sex, batch, and library metrics; missing covariates are not synthesized.

## Notes

Wave 2 establishes reproducible CN-expression association outputs for both normalization tracks.
Wave 3 extends RNA integration with EBV variant→expression associations and QC aggregation.

If `--wgs_only true` is used in `RNA_INTEGRATION`, RNA association steps in this contract are skipped by design.
