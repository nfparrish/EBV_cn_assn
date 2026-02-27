# HPC portability and deployment checklist (Phase 6 Wave 3)

## Overview

This document provides scheduler-portable execution guidance for the EBV multi-omics pipeline.

## Profile usage

The pipeline ships with baseline profiles in [pipeline/conf/base.config](pipeline/conf/base.config):
- `standard` (local executor)
- `slurm`
- `pbs`
- `dcc` (Duke DCC-oriented Slurm defaults)

Run examples:

### Local

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile standard`

### Slurm

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile slurm`

### PBS

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile pbs`

### Duke DCC (recommended)

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile dcc --slurm_account <your_account> --slurm_partition common --hpc_work_dir /work/<netid>/gsd/work --hpc_singularity_cache /work/<netid>/gsd/singularity-cache --outdir /work/<netid>/gsd/results`

## Deployment checklist

- Confirm Nextflow version and Java runtime on target cluster.
- Confirm container runtime policy (Singularity/Apptainer where required).
- Validate data paths and shared filesystem visibility across compute nodes.
- Run pilot manifest first (`pipeline/tests/fixtures/pilot_1000g_manifest.csv`).
- Capture report reproducibility outputs for each release run.

## DCC-specific guidance

- Prefer `/work` for `--outdir`, Nextflow `work` directory, and container cache to avoid filling `/hpc/home` (25 GB quota).
- Use `common` for general CPU runs when available to all accounts; use lab partitions when required by your allocation.
- Set `--slurm_account` explicitly for consistent charging and scheduling behavior.
- Start with workflow-specific runs (`intake`, `cn`, `variant`, `assoc`, `wgs_only`, `rna`) before full `report` on large cohorts.
- Keep walltimes under partition limits (`common` up to 90d default; `gpu-common` up to 2d; `scavenger-gpu` up to 7d).

## DCC examples

Pilot on CPU partition:

`nextflow run pipeline/main.nf --workflow report --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/host_variants.tsv -profile dcc --slurm_account <your_account> --slurm_partition common --hpc_work_dir /work/<netid>/gsd/work --hpc_singularity_cache /work/<netid>/gsd/singularity-cache --outdir /work/<netid>/gsd/results/pilot_report`

RNA FASTQ recount on CPU partition (RNA-enabled mode):

`nextflow run pipeline/main.nf --workflow rna --manifest <manifest.csv> --rna_input_mode fastq --rna_reference_fasta <combined_host_ebv.fa> --rna_reference_gtf <combined_host_ebv.gtf> -profile dcc --slurm_account <your_account> --slurm_partition common --hpc_work_dir /work/<netid>/gsd/work --hpc_singularity_cache /work/<netid>/gsd/singularity-cache --outdir /work/<netid>/gsd/results/rna_fastq`

Use a different partition:

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile dcc --slurm_account <your_account> --slurm_partition <allowed_partition> --slurm_qos <qos_if_required> --hpc_work_dir /work/<netid>/gsd/work --hpc_singularity_cache /work/<netid>/gsd/singularity-cache --outdir /work/<netid>/gsd/results/release`

## DCC sbatch wrapper (10 WGS test samples)

Use [pipeline/scripts/dcc_wgs10_testset.sbatch](pipeline/scripts/dcc_wgs10_testset.sbatch) to run an end-to-end test that:
- discovers paired WGS FASTQ runs from ENA study `PRJNA41139`
- downloads 10 samples to `/work`
- builds a valid manifest and host variant file
- launches `--workflow report` with the `dcc` profile

Submit example:

`sbatch --account <your_account> --partition common pipeline/scripts/dcc_wgs10_testset.sbatch`

Optional overrides at submit time:

`sbatch --account <your_account> --partition common --export=ALL,TESTSET_SIZE=10,STUDY_ACCESSION=PRJNA41139,WORK_BASE=/work/<netid>/gsd_wgs10,QOS=<optional_qos> pipeline/scripts/dcc_wgs10_testset.sbatch`

## Resource tuning

Edit process labels in [pipeline/conf/base.config](pipeline/conf/base.config):
- `manifest_validation`
- `cn_compute`
- `variant_compute`
- `assoc_compute`
- `rna_compute`

Tune `cpus`, `memory`, and `time` per site policy.

## Recommended release command

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> --outdir <release_outdir> -profile <standard|slurm|pbs|dcc>`
