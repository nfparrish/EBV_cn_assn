---
phase: "05"
name: "RNA-seq Integration"
created: 2026-02-26
---

# Phase 5: RNA-seq Integration — Context

## Decisions

- Phase 5 consumes completed Phase 2 CN outputs, Phase 3 variant outputs, and Phase 4 association infrastructure.
- Default RNA-seq source context is **MAGE-linked public RNA-seq**.
- Expression normalization baseline uses **two parallel outputs**:
  - `log2(TPM+1)`
  - quantile-normalized expression
- Baseline CN→expression association model is **linear model on normalized expression**.
- Baseline covariates include:
  - PC1-PC5
  - sex
  - batch
  - library metrics (e.g., depth/alignment-rate style technical covariates)
- Implementation remains Nextflow DSL2 + profile-driven HPC portability.

## Discretion Areas

- Exact ingestion schema for MAGE-linked manifests and optional matrix imports.
- How dual-normalization outputs are represented (long/wide) as long as both are emitted and documented.
- Default library metric set and fallback behavior when metrics are missing.
- Ordering and chunking strategy for CN-expression and variant-expression association execution.

## Deferred Ideas

- Mixed-model RNA analyses and random-effects extensions (v2 methods).
- Final manuscript/dashboard packaging (Phase 6).
- New cohort expansion beyond baseline source until Phase 5 baseline is stable.

## Inputs Required for Planning

- RNA-seq sample metadata and expression source tables (MAGE-linked public references).
- Phase 2 CN artifacts (`ebv_cn_normalized.tsv`, confidence/QC outputs).
- Phase 3 filtered EBV/host variant outputs.
- Covariate and technical metrics aligned to RNA sample IDs.

## Phase Boundary

Phase 5 delivers RNA-seq ingestion/normalization and associations for:
1. EBV CN → EBV gene expression
2. EBV variants → EBV gene expression

It does not finalize publication packaging.
