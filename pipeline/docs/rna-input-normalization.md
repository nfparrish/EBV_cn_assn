# RNA input and normalization contract (Phase 5 Wave 1)

## Overview

Phase 5 Wave 1 introduces RNA-seq ingestion and dual normalization outputs.

Run with:

`nextflow run pipeline/main.nf --workflow wgs_only --manifest <manifest.csv> -profile standard`

Optional inputs:

- `--rna_expression_input <expression.tsv>`
- `--rna_metadata_input <rna_metadata.tsv>`
- `--rna_require_ebv_genes <true|false>` (default: `true`)
- `--wgs_only <true|false>` (default: `false`)

WGS-only mode:
- Set `--workflow wgs_only` to run EBV depth/CN + EBV variant calling only.
- In this mode, RNA ingest/normalization/association modules are skipped.
- EBV outputs remain available via workflow emits (`ebv_depth_files`, `ebv_variant_tsv`, CN outputs).

RNA-enabled mode:
- Use `--workflow rna` when RNA ingest/normalization/association should execute.

Input mode:
- `--rna_input_mode auto|matrix|fastq` (default: `auto`)
	- `matrix`: use provided/derived matrix ingestion path
	- `fastq`: force raw FASTQ alignment + recount path
	- `auto`: use `fastq` when no `--rna_expression_input` is set and both reference files are provided; otherwise `matrix`

FASTQ recount mode requires:
- `--rna_reference_fasta <combined_host_ebv.fa>`
- `--rna_reference_gtf <combined_host_ebv.gtf>`
- `minimap2`, `samtools`, and `featureCounts` available on `PATH`

No deterministic fallback artifacts are generated.
Matrix mode is empirical-only: `--rna_expression_input` is required and synthetic expression/metadata values are never generated.

## Ingestion outputs

- `rna_expression_canonical.tsv`
- `rna_sample_map.tsv`
- `rna_ingest_report.json`
- `rna_fastq_count_report.json` (only in FASTQ recount mode)

When `--wgs_only true`, RNA ingestion outputs are not produced.

By default, ingest enforces EBV expression presence. If no EBV-like gene rows are detected in the canonical matrix, the ingest step fails unless `--rna_require_ebv_genes false` is set.

### `rna_expression_canonical.tsv`

Expected shape:

- rows: genes (`gene_id`)
- columns: RNA sample IDs (aligned to manifest sample IDs)

### `rna_sample_map.tsv`

Fields:

- `sample_id`, `dataset`, `population`, `sex`, `sequence_lab`
- `rna_batch`, `library_size`, `mapping_rate`
- `rna_source`

## Normalization outputs

- `rna_expression_log2tpm.tsv` (`log2(TPM+1)`)
- `rna_expression_quantile.tsv` (quantile normalization)
- `rna_normalization_report.json`

Both normalized outputs preserve canonical gene/sample dimensions.

When `--wgs_only true`, RNA normalization and RNA association outputs are skipped.

## Input expectations

### Expression input (optional)

Tab-delimited table with:

- first column: `gene_id`
- remaining columns: sample IDs matching manifest `sample_id` values

### Metadata input (optional)

Tab-delimited table with at least:

- `sample_id`

Optional metadata fields used when present:

- `rna_batch`, `library_size`, `mapping_rate`, `sex`, `sequence_lab`

No random or synthetic fallback values are generated for missing metadata fields.

## EBV coverage and ingest report fields

`rna_ingest_report.json` includes:

- `ebv_gene_count`
- `non_ebv_gene_count`
- `ebv_gene_required`

## Notes

- Wave 1 focuses on reproducible ingestion + normalization contracts.
- Association execution against RNA outputs is introduced in later Phase 5 waves.
