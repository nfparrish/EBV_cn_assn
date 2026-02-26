# RNA FASTQ recount contract

This mode enables direct ingestion of raw RNA-seq FASTQ files, alignment to a combined host+EBV reference, and generation of a fresh gene count matrix.

## Purpose
- Intake raw RNA-seq reads (`read1_path`/`read2_path` from manifest)
- Map to reference containing both host and EBV contigs
- Generate gene-level count table for downstream normalization and association

## Required parameters
- `--workflow rna`
- `--rna_input_mode fastq` (or `auto` with references and no `rna_expression_input`)
- `--rna_reference_fasta <combined_host_ebv.fa>`
- `--rna_reference_gtf <combined_host_ebv.gtf>`

Manifest FASTQ columns used in this mode:
- preferred: `rna_read1_path`, `rna_read2_path`
- fallback: `read1_path`, `read2_path`

## Required tools
- `minimap2`
- `samtools`

Counter backends:
- `featureCounts` (if available)
- built-in fallback counter (BAM + GTF overlap counting) when `featureCounts` is unavailable

Use `--rna_counter auto` (default) to auto-select best available backend.

## Combined reference validation
The recount step validates that references include both host and EBV content:
- FASTA must contain at least one EBV contig (e.g. `NC_007605`) and at least one non-EBV contig.
- GTF must include EBV-annotated genes plus host genes.

## Outputs
- `rna_expression_canonical.tsv`
- `rna_sample_map.tsv`
- `rna_ingest_report.json`
- `rna_fastq_count_report.json`

### `rna_fastq_count_report.json`
Includes:
- toolchain used (`aligner`, `counter`)
- reference paths
- sample/gene counts
- per-sample `library_size` and `mapping_rate`

## Notes
- FASTQ recount expects local FASTQ paths in the manifest.
- Existing matrix-input mode remains supported and unchanged.
