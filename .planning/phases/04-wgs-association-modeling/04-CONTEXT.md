---
phase: "04"
name: "WGS Association Modeling"
created: 2026-02-26
---

# Phase 4: WGS Association Modeling — Context

## Decisions

- Phase 4 consumes Phase 2 CN outputs and Phase 3 filtered/provenance variant outputs.
- Baseline model family is **generalized linear model (GLM)**.
- Outcome handling baseline is **z-scored EBV CN** before model fitting.
- Multiple-testing correction baseline is **FDR (Benjamini-Hochberg)**.
- Covariates should use **auto fallback**:
  - use separate `sequence_lab` and `batch` if both present;
  - otherwise derive/use a combined `lab_batch` style factor.
- Required covariate intent remains: PCs (PC1-PC5), sex, and sequencing/lab batch structure.
- Implementation remains Nextflow DSL2 + profile-driven HPC portability.

## Discretion Areas

- Exact GLM family/link defaults (e.g., gaussian/identity baseline with explicit overrides).
- Variant-wise model chunking and parallelization strategy.
- Diagnostics thresholds and plotting detail depth.
- Output table partitioning for EBV and host variant tracks.

## Deferred Ideas

- Mixed-effects and relatedness-aware models (v2 method expansion).
- RNA-seq expression association workflows (Phase 5).
- Final publication/dashboard packaging (Phase 6).

## Inputs Required for Planning

- `ebv_cn_normalized.tsv` and confidence/QC artifacts (Phase 2).
- Filtered EBV and host variant outputs + provenance artifacts (Phase 3).
- Covariate table with PC1-PC5, sex, and sequencing metadata.

## Phase Boundary

Phase 4 delivers covariate-harmonized, variant-wise association statistics with FDR-adjusted p-values and model diagnostics. It does not include RNA-seq analyses.
