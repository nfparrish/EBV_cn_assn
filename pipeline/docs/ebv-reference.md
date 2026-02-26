# EBV reference preparation

## Overview

The EBV CN workflow supports two reference selection modes:

1. **Default accession mode** (`params.viral_reference_fasta` unset)
   - Uses `params.reference_accession` (default: `NC_007605`).
2. **User FASTA mode** (`params.viral_reference_fasta` set)
   - Uses the user-provided viral FASTA path and records that source.

`EBV_REFERENCE_PREP` emits `ebv_reference_manifest.json` as the deterministic contract consumed by downstream CN modules.

## Output contract

`ebv_reference_manifest.json`

- `reference_accession`: configured accession identifier.
- `source_type`: one of `default_accession` or `user_fasta`.
- `source`: accession ID or FASTA path depending on mode.

## Reproducibility notes

- Reference choice is encoded explicitly in the manifest and propagated to per-sample and cohort outputs.
- All downstream outputs include `reference_accession` for traceability.
