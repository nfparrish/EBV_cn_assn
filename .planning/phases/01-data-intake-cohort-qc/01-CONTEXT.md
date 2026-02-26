---
phase: "01"
name: "Data Intake & Cohort QC"
created: 2026-02-26
---

# Phase 1: Data Intake & Cohort QC — Context

## Decisions

- Use Nextflow DSL2 with nf-core style module boundaries for all phase-1 implementation.
- Treat 1000 Genomes Phase 3 as the required initial WGS source for v1.
- Pipeline must be scheduler-agnostic for HPC portability (profile-driven, not cluster hard-coded).
- Phase-1 outputs must include per-sample QC and cohort-level summary artifacts that feed downstream CN/association phases.
- Input handling must be resumable and support restart without reprocessing successful samples.

## Discretion Areas

- Exact manifest schema details (column names and optional metadata) as long as required sample identifiers and file locations are present.
- Choice of concrete QC metric set and threshold strategy for initial implementation, provided outputs are reproducible and interpretable.
- Directory layout for pipeline code (`main.nf`, `modules/`, `conf/`, `assets/`) as long as it remains nf-core compatible.
- Container strategy (Docker/Singularity/Conda) for local vs HPC execution profiles.

## Deferred Ideas

- Large-scale association model execution (covered in Phase 4).
- RNA-seq ingestion and expression association workflow (covered in Phase 5).
- Interactive dashboard finalization and publication-ready reporting polish (covered in Phase 6).
- Mixed-model and rare-variant methods (v2 scope).

## Inputs Required for Planning

- WGS source baseline: `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/`
- Required covariate awareness downstream: PC1-PC5, sex, sequencing lab/batch.

## Phase Boundary

Phase 1 delivers only reliable data intake and quality-control foundations. It does not compute EBV copy number yet, but it must produce validated sample metadata and QC artifacts that make CN computation and statistical modeling straightforward in subsequent phases.

## Pre-Implementation Checklist

- [ ] Confirm pilot manifest subset for smoke tests (coverage/diversity representative).
- [ ] Confirm container strategy per environment (Docker local, Singularity/Apptainer on HPC).
- [ ] Confirm initial QC policy (warning-first thresholds, hard filters deferred).
- [ ] Confirm scheduler profiles to support first (at minimum base + Slurm; add PBS profile if required now).
- [ ] Confirm output contract for Phase 1 artifacts (sample status table, per-sample QC files, cohort summary report).
