---
phase: 03-variant-calling-layer
plan: "02"
subsystem: host-variant-ingest-harmonize
tags:
  - nextflow
  - host-variants
  - harmonization
  - association-input
provides:
  - Host variant ingestion from public variant resource inputs
  - Harmonized host variant table aligned to manifest cohort metadata
affects:
  - 03-01
tech-stack:
  added:
    - Nextflow DSL2 modules: HOST_VARIANT_INGEST, HOST_VARIANT_HARMONIZE
    - Python scripts for host input staging and cohort harmonization
  patterns:
    - Staged-input validation + explicit sample overlap/mismatch reporting
key-files:
  created:
    - pipeline/modules/local/host_variant_ingest/main.nf
    - pipeline/modules/local/host_variant_harmonize/main.nf
    - pipeline/scripts/stage_host_variant_input.py
    - pipeline/scripts/harmonize_host_variants.py
    - pipeline/docs/host-variant-inputs.md
    - pipeline/tests/fixtures/pilot_host_variants.tsv
  modified:
    - pipeline/workflows/ebv_variant_layer.nf
    - pipeline/conf/base.config
key-decisions:
  - decision: Require `--host_variant_input` for `--workflow variant` runs in Wave 2+.
    rationale: Enforces VAR-02 readiness and prevents silent execution without host variant pathway coverage.
  - decision: Emit explicit mismatch artifacts for input-vs-manifest sample IDs.
    rationale: Keeps host pathway auditable and directly useful for downstream association matrix QC.
patterns-established:
  - Host variant path now runs in parallel with EBV variant path and emits harmonized association-ready outputs.
duration: 45min
completed: 2026-02-26
---

# Phase 3: Variant Calling Layer Summary (03-02)

Implemented host variant ingestion and harmonization pathway for Wave 2.

## Accomplishments
- Added `HOST_VARIANT_INGEST` and `HOST_VARIANT_HARMONIZE` modules and wired them into `EBV_VARIANT_LAYER`.
- Added host variant config params (`host_variant_input`, `host_variant_format`) and required host input guard for `variant` workflow runs.
- Added host variant staging/harmonization scripts with sample-overlap and mismatch reporting.
- Added host variant operator documentation and pilot host fixture.

## Output contracts added
- `host_variant_staged.tsv`
- `host_variant_sample_map.tsv`
- `host_variant_ingest_report.json`
- `host_variant_harmonized.tsv`
- `host_variant_harmonize_report.json`
- `host_variant_sample_mismatch.tsv`

## Validation Notes
- Successful pilot run:
  `nextflow run pipeline/main.nf --workflow variant --manifest pipeline/tests/fixtures/pilot_1000g_manifest.csv --host_variant_input pipeline/tests/fixtures/pilot_host_variants.tsv --host_variant_format tsv --outdir pipeline/tests/pilot_variant_run_03_02 -profile standard -resume`
- Host ingest and harmonization completed successfully alongside EBV variant modules.

## Next Readiness
- Phase 3 Wave 3 (`03-03`) can now implement balanced filtering + provenance reporting across both EBV and host variant pathways.
