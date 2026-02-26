# EBV variant calling (Phase 3 Wave 1)

## Overview

Phase 3 Wave 1 implements EBV variant generation and cohort matrix aggregation with a Mutect2 baseline strategy and dual matrix outputs.

Workflow path:

1. `EBV_REFERENCE_PREP` (reference manifest contract)
2. `EBV_VARIANT_CALL` (per-sample EBV variants)
3. `EBV_VARIANT_MATRIX` (cohort AF and binary matrices)

Run with:

`nextflow run pipeline/main.nf --workflow variant --manifest <manifest.csv> -profile standard`

## Per-sample call artifact

Each sample emits `<sample_id>.ebv_variants.tsv` with:

- sample metadata (`sample_id`, `dataset`, `population`, `sex`, `sequence_lab`)
- variant fields (`chrom`, `pos`, `ref`, `alt`, `af`, `dp`, `qual`, `filter`)
- provenance (`caller`, `reference_accession`)

Caller provenance currently records `mutect2_baseline` for downstream tracking.

## Cohort outputs

- `ebv_variant_catalog.tsv`: variant list with carrier counts and mean AF.
- `ebv_variant_af_matrix.tsv`: sample x variant AF matrix.
- `ebv_variant_binary_matrix.tsv`: sample x variant binary matrix using AF threshold.
- `ebv_variant_matrix_summary.json`: matrix dimensions and filter-level aggregate counts.

## Current defaults

- `params.variant_caller = mutect2`
- `params.variant_filter_profile = balanced`
- `params.variant_binary_af_threshold = 0.05`

## Notes

Wave 1 establishes the EBV variant output contract for downstream filtering/provenance work in `03-03` and host harmonization integration in `03-02`.
