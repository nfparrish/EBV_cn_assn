# EBV Copy Number Multi-Omics Pipeline

## What This Is

An HPC-portable analysis pipeline for large public datasets to estimate Epstein-Barr virus (EBV) copy number from whole-genome sequencing (WGS) data, then test associations between EBV copy number and both viral and human genetic variation. The project also extends to public RNA-seq datasets to evaluate associations between EBV copy number/EBV variants and EBV gene expression.

## Core Value

Deliver a reproducible, scheduler-agnostic pipeline that can process thousands of public samples and produce statistically robust EBV copy number association results.

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] Process WGS cohorts at scale and compute per-sample EBV copy number.
- [ ] Call EBV and host variants and test associations with EBV copy number using covariate-adjusted models.
- [ ] Integrate RNA-seq analysis for EBV expression associations and output publication-ready artifacts.

### Out of Scope

- Clinical/diagnostic deployment — this project is research-focused on public datasets.
- Interactive web app for non-technical users — v1 prioritizes reproducible pipeline execution and reports.

## Context

- Primary WGS source for v1: 1000 Genomes Phase 3 (`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/`).
- RNA-seq reference requested for v1 context: `https://github.com/mccoy-lab/MAGE`.
- Preferred pipeline framework: Nextflow with nf-core style conventions.
- Statistical baseline: linear models with covariates (first 5 genetic PCs, sex, sequencing lab/batch).
- Desired outputs include QC reports, per-sample EBV CN calls, association result tables, manuscript-ready figures, and an interactive dashboard.

## Constraints

- **Tech stack**: Nextflow + nf-core style workflow — required for reproducibility and portability.
- **Compute environment**: Scheduler-agnostic HPC deployment — must run with profile-based support for major schedulers.
- **Scale**: Thousands of public WGS samples — requires parallelized, restartable workflow design.
- **Data scope**: Publicly available datasets only for v1 — simplifies legal and access constraints.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Use Nextflow (nf-core style) | Standardized, portable HPC execution and reproducibility | — Pending |
| Start with linear association models | Fast baseline with interpretable covariate adjustment | — Pending |
| Target scheduler-agnostic profiles | Pipeline should transfer cleanly across HPC environments | — Pending |

---
*Last updated: 2026-02-26 after project initialization*
