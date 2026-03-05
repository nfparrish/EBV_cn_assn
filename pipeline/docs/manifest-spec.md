# WGS Manifest Specification (Phase 1)

This manifest defines cohort input for Phase 1 data intake and QC.

## Required Columns

- `sample_id` — unique sample identifier; `[A-Za-z0-9._-]+`
- `dataset` — `1000G_PHASE3` (preferred for v1) or `OTHER_PUBLIC`
- `population` — cohort/population label
- `sex` — `M`, `F`, or `U`
- `sequence_lab` — sequencing site/batch origin label

## Sequencing Input Columns

Each row must include at least one of:

1. paired FASTQ: `read1_path` + `read2_path`
2. aligned BAM: `bam_path`
3. aligned CRAM: `cram_path`

Optional:

- `md5` — checksum string for source file validation
- `notes` — free-text provenance annotations
- `rna_read1_path`, `rna_read2_path` — optional RNA-specific FASTQ paths used by RNA recount workflows; when present, these are preferred over `read1_path`/`read2_path` for RNA ingestion.
- `family_id` — pedigree family identifier (matches `.ped` file Family ID)
- `father_id` — sample_id of father (`0` if unknown)
- `mother_id` — sample_id of mother (`0` if unknown)
- `relationship` — `father`, `mother`, `child`, or `unrelated`

These pedigree columns enable trio identification for parent-offspring heritability estimation without a separate sidecar mapping.

## Behavior

- Validation is strict and fail-fast.
- Duplicate `sample_id` values are rejected.
- Invalid rows are reported in `manifest_validation_report.json`.
- On success, rows are normalized and written to `validated_manifest.csv`.

## Example

```csv
sample_id,dataset,population,sex,sequence_lab,read1_path,read2_path,bam_path,md5
HG00096,1000G_PHASE3,GBR,M,LAB_A,/data/HG00096_R1.fastq.gz,/data/HG00096_R2.fastq.gz,,0123456789abcdef0123456789abcdef
HG00097,1000G_PHASE3,GBR,F,LAB_A,,,/data/HG00097.bam,
```

For aligned-only rows, provide `bam_path` or `cram_path`.
