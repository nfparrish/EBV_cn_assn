---
phase: 04-wgs-association-modeling
plan: "02"
subsystem: association-execution
tags:
  - nextflow
  - association
  - glm
  - fdr
provides:
  - EBV and host variant-wise association result tables
  - Unified inferential schema with FDR-adjusted p-values
affects:
  - 04-01
tech-stack:
  added:
    - Nextflow DSL2 modules: ASSOC_RUN_EBV, ASSOC_RUN_HOST
    - Shared Python association runner with covariate-adjusted residualized regression and BH-FDR
  patterns:
    - Unified EBV/host result schema and pathway-specific reports
key-files:
  created:
    - pipeline/modules/local/assoc_run_ebv/main.nf
    - pipeline/modules/local/assoc_run_host/main.nf
    - pipeline/scripts/run_assoc_variants.py
    - pipeline/docs/association-output-contract.md
  modified:
    - pipeline/workflows/wgs_association.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Reuse Phase 3 filtered variant outputs directly in association workflow.
    rationale: Ensures Wave 2 statistics align with existing VAR-03 filtering/provenance policy.
  - decision: Emit strict JSON reports (`null` for non-finite minima) and standardized result columns.
    rationale: Keeps downstream parsing stable and machine-readable.
patterns-established:
  - `--workflow assoc` now runs CN prep, variant layer, and EBV/host association execution in one reproducible chain.
duration: 50min
completed: 2026-02-26
---

# Phase 4: WGS Association Modeling Summary (04-02)

Implemented Wave 2 association execution for EBV and host variant pathways.

## Accomplishments
- Added `ASSOC_RUN_EBV` and `ASSOC_RUN_HOST` modules.
- Added shared association runner script emitting effect size, standard error, z-score, p-value, and FDR-adjusted p-value.
- Wired `WGS_ASSOCIATION` to consume Phase 3 filtered variant artifacts (`ebv_variant_filtered`, `host_variant_filtered`).
- Added output contract documentation for association result artifacts and reports.
- Added `association_min_samples` configuration parameter.

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow assoc --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants.tsv --host_variant_format tsv --outdir pipeline/tests/pilot_assoc_run_04_02b -profile standard`
- Both `ASSOC_RUN_EBV` and `ASSOC_RUN_HOST` completed successfully.
- Result tables include required inferential columns and FDR fields.

## Next Readiness
- Wave 3 (`04-03`) can now add diagnostics and association QC summary reporting on top of generated EBV/host results.
