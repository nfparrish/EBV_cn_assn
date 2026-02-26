---
status: complete
phase: 01-data-intake-cohort-qc
source: 01-01-SUMMARY.md, 01-02-SUMMARY.md, 01-03-SUMMARY.md
started: 2026-02-26T12:30:00Z
updated: 2026-02-26T12:49:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Manifest Validation Gate
expected: Invalid manifest rows fail fast and valid rows normalize into validated manifest output.
result: pass

### 2. Resumable Status Registry
expected: Sample status summary and registry are emitted with processed/skipped counts for rerun control.
result: pass

### 3. Per-Sample QC Artifacts
expected: One QC artifact is generated per sample with metrics and warning fields.
result: pass

### 4. Cohort QC Summary Bundle
expected: Cohort-level machine-readable outputs and markdown report are generated.
result: pass

### 5. Local Pilot Execution
expected: Phase-1 pipeline runs locally on a small two-sample 1000G-based manifest and produces final artifacts.
result: pass

## Summary

total: 5
passed: 5
issues: 0
pending: 0
skipped: 0

## Notes

- Pilot manifest: `pipeline/tests/fixtures/pilot_1000g_manifest.csv`
- Pilot artifacts: `pipeline/tests/pilot_run/artifacts/`
- Observed pilot summary (`qc_cohort_summary.json`):
  - sample_count: 2
  - processed: 2
  - skipped: 0
  - q30_rate_mean: 0.878
  - mapping_rate_mean: 0.9578
  - mean_coverage_mean: 29.81
  - low_q30 warnings: 1

## Gaps

None.
