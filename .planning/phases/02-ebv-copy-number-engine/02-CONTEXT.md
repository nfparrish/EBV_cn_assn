---
phase: "02"
name: "EBV Copy Number Engine"
created: 2026-02-26
---

# Phase 2: EBV Copy Number Engine — Context

## Decisions

- Phase 2 must consume Phase 1 intake outputs (validated sample metadata and QC context) as upstream inputs.
- EBV copy number is the primary outcome and must be produced per sample (CN-01).
- Cohort output must include normalized EBV CN table with metadata fields required for association modeling (CN-02).
- Confidence and quality annotations are required for each CN estimate (CN-03).
- Implementation remains Nextflow DSL2 + profile-driven HPC portability.

## Discretion Areas

- Exact EBV reference asset packaging and indexing strategy.
- Specific CN estimation method details (read-depth normalization formulation) as long as outputs are reproducible and auditable.
- Confidence scoring approach and threshold defaults, provided they are documented and emitted with artifacts.
- Internal file layout for CN modules/workflows as long as it fits existing `pipeline/` conventions.

## Deferred Ideas

- Host/EBV variant association analysis (Phase 4).
- RNA-seq expression association analysis (Phase 5).
- Final reporting/dashboard polish (Phase 6).
- Mixed-model / rare-variant methods (v2).

## Inputs Required for Planning

- Upstream intake/QC artifacts from Phase 1.
- 1000 Genomes pilot subset for method validation and quality checks.

## Phase Boundary

Phase 2 delivers EBV copy number computation, normalization-ready cohort outputs, and confidence/QC annotations. It does not run association models yet.

