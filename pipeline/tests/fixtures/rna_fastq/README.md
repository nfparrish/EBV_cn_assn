# RNA FASTQ fixture

This fixture provides a minimal local dataset for testing FASTQ-based RNA recount mode.

## Files
- `pilot_rna_fastq_manifest.csv` — manifest with local FASTQ paths.
- `combined_host_ebv.fa` — tiny combined host+EBV reference FASTA.
- `combined_host_ebv.gtf` — tiny combined host+EBV annotation GTF.
- `HG00096_R1.fastq`, `HG00096_R2.fastq`
- `HG00097_R1.fastq`, `HG00097_R2.fastq`

## Example run

```bash
nextflow run pipeline/main.nf \
  --workflow rna \
  --manifest pipeline/tests/fixtures/rna_fastq/pilot_rna_fastq_manifest.csv \
  --rna_input_mode fastq \
  --rna_reference_fasta pipeline/tests/fixtures/rna_fastq/combined_host_ebv.fa \
  --rna_reference_gtf pipeline/tests/fixtures/rna_fastq/combined_host_ebv.gtf \
  --rna_assoc_min_samples 2 \
  --outdir pipeline/tests/pilot_rna_fastq_run \
  -profile standard
```

Prerequisites on PATH: `minimap2`, `samtools`, `featureCounts`.

If `featureCounts` is not installed, keep `--rna_counter auto` (default) and the workflow will use the built-in counting fallback.
