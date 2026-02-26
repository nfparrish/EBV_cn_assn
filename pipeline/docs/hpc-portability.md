# HPC portability and deployment checklist (Phase 6 Wave 3)

## Overview

This document provides scheduler-portable execution guidance for the EBV multi-omics pipeline.

## Profile usage

The pipeline ships with baseline profiles in [pipeline/conf/base.config](pipeline/conf/base.config):
- `standard` (local executor)
- `slurm`
- `pbs`

Run examples:

### Local

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile standard`

### Slurm

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile slurm`

### PBS

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> -profile pbs`

## Deployment checklist

- Confirm Nextflow version and Java runtime on target cluster.
- Confirm container runtime policy (Singularity/Apptainer where required).
- Validate data paths and shared filesystem visibility across compute nodes.
- Run pilot manifest first (`pipeline/tests/fixtures/pilot_1000g_manifest.csv`).
- Capture report reproducibility outputs for each release run.

## Resource tuning

Edit process labels in [pipeline/conf/base.config](pipeline/conf/base.config):
- `manifest_validation`
- `cn_compute`
- `variant_compute`
- `assoc_compute`
- `rna_compute`

Tune `cpus`, `memory`, and `time` per site policy.

## Recommended release command

`nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host.tsv> --outdir <release_outdir> -profile <standard|slurm|pbs>`
