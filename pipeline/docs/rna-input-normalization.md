# RNA input and normalization contract (Phase 5 Wave 1)

## Overview

Phase 5 Wave 1 introduces RNA-seq ingestion and dual normalization outputs.

Run with:

`nextflow run pipeline/main.nf --workflow rna --manifest <manifest.csv> -profile standard`

Optional inputs:

- `--rna_expression_input <expression.tsv>`
- `--rna_metadata_input <rna_metadata.tsv>`
- `--rna_require_ebv_genes <true|false>` (default: `true`)

Input mode:
- `--rna_input_mode auto|matrix|fastq` (default: `auto`)
	- `matrix`: use provided/derived matrix ingestion path
	- `fastq`: force raw FASTQ alignment + recount path
	- `auto`: use `fastq` when no `--rna_expression_input` is set and both reference files are provided; otherwise `matrix`

FASTQ recount mode requires:
- `--rna_reference_fasta <combined_host_ebv.fa>`
- `--rna_reference_gtf <combined_host_ebv.gtf>`
- `minimap2`, `samtools`, and `featureCounts` available on `PATH`

If optional inputs are omitted, deterministic fallback artifacts are generated for pilot execution.

## Ingestion outputs

- `rna_expression_canonical.tsv`
- `rna_sample_map.tsv`
- `rna_ingest_report.json`
- `rna_fastq_count_report.json` (only in FASTQ recount mode)

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

## EBV coverage and ingest report fields

`rna_ingest_report.json` includes:

- `ebv_gene_count`
- `non_ebv_gene_count`
- `ebv_gene_required`

## Notes

- Wave 1 focuses on reproducible ingestion + normalization contracts.
- Association execution against RNA outputs is introduced in later Phase 5 waves.
