# RNA Association + QC Contract

This document defines Phase 5 Wave 3 outputs for EBV variant→expression association and integrated RNA association QC.

## Variant→Expression association outputs

Process: `RNA_ASSOC_VARIANT`

### `rna_assoc_variant_results.tsv`
Tab-delimited schema:
- `normalization` — normalization label (`log2_tpm_plus_1` or `quantile_normalized`)
- `variant_id` — EBV variant key (`chrom:pos:ref>alt`)
- `gene_id` — EBV gene/feature identifier
- `n_samples` — overlap sample count used for fit
- `n_nonzero_predictor` — number of samples with non-zero variant predictor value
- `effect_size` — linear model coefficient for variant predictor
- `std_error` — standard error of effect estimate
- `t_stat` — Wald-style test statistic approximation
- `p_value` — nominal two-sided p-value
- `p_value_fdr` — BH-adjusted p-value across all variant-gene tests in the output
- `model_family` — model family label (`gaussian` baseline)
- `model_link` — link label (`identity` baseline)
- `covariate_terms` — semicolon-delimited covariate term list used in fit

### `rna_assoc_variant_report.json`
Summary contract:
- `normalization`
- `variant_input_row_count`
- `variant_count_modeled`
- `test_count`
- `gene_count`
- `sample_overlap_count`
- `min_samples_required`
- `covariate_terms`
- `filters.min_af`
- `filters.min_dp`
- `filters.min_qual`
- `min_p_value`

## RNA integration QC outputs

Process: `RNA_ASSOC_QC`

### `rna_assoc_qc_summary.json`
- `coverage.sample_overlap` per association layer
- `coverage.gene_counts` per association layer
- `coverage.variant_count_modeled` for variant layers
- `significance` block per layer:
  - `test_count`
  - `finite_p_count`
  - `finite_p_rate`
  - `modeled_count`
  - `modeled_rate`
  - `nominal_hit_count`
  - `fdr_hit_count`
- `model_health` high-level flags:
  - `all_layers_have_finite_p`
  - `low_coverage_any_layer`
  - `all_layers_nominally_null`

### `rna_assoc_qc_report.md`
Human-readable markdown summary with coverage, significance, and model-health flags.

### `rna_assoc_qc_top_hits.tsv`
Top-ranked hits across CN and variant layers:
- `layer`
- `feature_id`
- `gene_id`
- `effect_size`
- `p_value`
- `p_value_fdr`
- `n_samples`
