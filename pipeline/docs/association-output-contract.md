# Association output contract (Phase 4 Wave 2)

## Overview

Phase 4 Wave 2 adds variant-wise association execution for both EBV and host pathways using a shared output schema and FDR-adjusted p-values.

Outputs are produced by:

- `ASSOC_RUN_EBV`
- `ASSOC_RUN_HOST`

## Required input artifacts

- `association_input.tsv` from Wave 1
- `ebv_variants_filtered.tsv` from Phase 3 filtering
- `host_variant_filtered.tsv` from Phase 3 filtering

## Result tables

### `assoc_ebv_results.tsv`
### `assoc_host_results.tsv`

Shared columns:

- `variant_id`
- `pathway` (`ebv` or `host`)
- `n_samples`
- `n_nonzero_predictor`
- `effect_size`
- `std_error`
- `z_score`
- `p_value`
- `p_value_fdr`
- `model_family`
- `model_link`
- `covariate_mode`

## Report JSON files

- `assoc_ebv_report.json`
- `assoc_host_report.json`

Report fields include pathway-specific result counts, covariate mode/terms, minimum p-value, and minimum sample threshold.

## Notes

- Multiple testing correction uses Benjamini-Hochberg FDR.
- Wave 2 establishes a unified inferential schema for EBV and host variant outputs.
