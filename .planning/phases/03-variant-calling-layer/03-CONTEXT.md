---
phase: "03"
name: "Variant Calling Layer"
created: 2026-02-26
---

# Phase 3: Variant Calling Layer — Context

## Decisions

- Phase 3 must consume Phase 2 EBV CN outputs and confidence context as upstream inputs.
- EBV variant calling baseline for v1 is **GATK Mutect2**.
- Host variant pathway for v1 is **consume validated public host VCF/PLINK** (not de novo host calling).
- Filtering posture is **balanced** (neither ultra-conservative nor permissive).
- EBV cohort outputs must include **both**:
  - allele-fraction matrix
  - derived binary presence/absence matrix
- Filtering provenance and retained/excluded counts are required for reproducibility (VAR-03).
- Implementation remains Nextflow DSL2 + profile-driven HPC portability.

## Discretion Areas

- Exact Mutect2 preprocessing granularity and calling windows.
- Specific balanced threshold defaults (QUAL/DP/AF), as long as they are configurable and reported.
- Final EBV matrix schema details (wide vs long companion formats), as long as both AF and binary outputs are emitted.
- Host variant harmonization internals (VCF vs PLINK normalization path) as long as association-ready outputs are produced.

## Deferred Ideas

- Association modeling itself (Phase 4).
- RNA-seq integration (Phase 5).
- Dashboard/publication packaging (Phase 6).
- Multi-caller ensemble and rare-variant burden workflows (v2+).

## Inputs Required for Planning

- Phase 2 outputs:
  - `ebv_cn_normalized.tsv`
  - `ebv_cn_confidence.tsv`
  - `ebv_cn_qc_summary.json`
- Pilot manifest from Phase 1 for run validation.
- Public host variant resources (VCF/PLINK) aligned to cohort sample IDs.

## Phase Boundary

Phase 3 delivers analysis-ready EBV and host variant inputs with filtering/provenance reports. It does not run association models.
