# EBV CN output contract

## Scope

Phase 2 Wave 1 emits raw and normalized EBV copy-number tables based on two methods:

- `depth_ratio_cn`
- `read_count_ratio_cn`

Host baseline used in this implementation is autosomal depth (`cn_host_baseline=autosomal`).

## Per-sample raw artifact

Each sample produces `<sample_id>.ebv_cn_raw.json` with:

- Sample metadata: `sample_id`, `dataset`, `population`, `sex`, `sequence_lab`
- Reference metadata: `reference_accession`
- Method values:
  - `methods.depth_ratio_cn`
  - `methods.read_count_ratio_cn`
- Formula trace:
  - `formula_trace.depth_ratio`
  - `formula_trace.read_count_ratio`
- Supporting metrics:
  - `host_autosomal_depth`, `ebv_depth`, `total_reads`, `ebv_mapped_reads`

## Cohort outputs

### `ebv_cn_raw_cohort.tsv`

Columns:

- `sample_id`
- `dataset`
- `population`
- `sex`
- `sequence_lab`
- `reference_accession`
- `depth_ratio_cn`
- `read_count_ratio_cn`
- `host_autosomal_depth`
- `ebv_depth`
- `total_reads`
- `ebv_mapped_reads`

### `ebv_cn_normalized.tsv`

Columns:

- `sample_id`
- `dataset`
- `population`
- `sex`
- `sequence_lab`
- `reference_accession`
- `depth_ratio_cn`
- `read_count_ratio_cn`
- `depth_ratio_cn_norm`
- `read_count_ratio_cn_norm`

Normalization in Wave 1 is cohort-mean scaling per method.

### `ebv_cn_cohort_summary.json`

Fields:

- `sample_count`
- `methods.depth_ratio_cn_mean`
- `methods.read_count_ratio_cn_mean`
- `host_baseline`
- `notes`

## Downstream expectation

Wave 2 confidence and QC modules consume `ebv_cn_normalized.tsv` and cohort summary statistics.
