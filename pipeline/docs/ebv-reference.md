# EBV reference preparation

## Overview

The EBV CN workflow supports two reference selection modes:

1. **Bundled FASTA default mode** (`params.viral_reference_fasta` unset)
   - Uses the repository-bundled NCBI EBV FASTA at `pipeline/data/reference/ebv/NC_007605.1.fa`.
   - `reference_accession` remains recorded (default: `NC_007605`).
2. **User FASTA mode** (`params.viral_reference_fasta` set)
   - Uses the user-provided viral FASTA path and records that source.

`EBV_REFERENCE_PREP` emits `ebv_reference_manifest.json` as the deterministic contract consumed by downstream CN modules.

## Output contract

`ebv_reference_manifest.json`

- `reference_accession`: configured accession identifier.
- `source_type`: one of `user_fasta` or `default_accession` (legacy).
- `source`: FASTA path (bundled default or user-provided), or accession ID in legacy mode.

## Reproducibility notes

- Reference choice is encoded explicitly in the manifest and propagated to per-sample and cohort outputs.
- All downstream outputs include `reference_accession` for traceability.
