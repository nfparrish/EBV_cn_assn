# Host variant inputs (Phase 3 Wave 2)

## Overview

Phase 3 Wave 2 adds a host variant ingestion and harmonization path that consumes validated public host variant resources and aligns them to cohort manifest sample IDs.

Run with:

`nextflow run pipeline/main.nf --workflow variant --manifest <manifest.csv> --host_variant_input <host_variant.tsv> --host_variant_format tsv -profile standard`

## Required parameters

- `--host_variant_input`: path to host variant input file.
- `--host_variant_format`: one of `tsv`, `vcf`, `plink`.

Current implementation expects a tabular export with required columns and records the declared format for provenance.

## TSV required columns

- `sample_id`
- `variant_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gt`
- `dosage`

## Module outputs

### `HOST_VARIANT_INGEST`

- `host_variant_staged.tsv`
- `host_variant_sample_map.tsv`
- `host_variant_ingest_report.json`

### `HOST_VARIANT_HARMONIZE`

- `host_variant_harmonized.tsv`
- `host_variant_harmonize_report.json`
- `host_variant_sample_mismatch.tsv`

## Harmonized schema

`host_variant_harmonized.tsv` includes:

- sample metadata from manifest (`dataset`, `population`, `sex`, `sequence_lab`)
- host variant fields (`variant_id`, `chrom`, `pos`, `ref`, `alt`, `gt`, `dosage`)

## Notes

- This Wave 2 path consumes host variants; it does not perform de novo host calling.
- Mismatch reporting is emitted for audit and downstream association readiness checks.
