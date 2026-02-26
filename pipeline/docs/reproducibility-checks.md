# Reporting reproducibility checks (Phase 6 Wave 3)

## Overview

Wave 3 adds reproducibility verification artifacts to `--workflow report` runs.

## Outputs

- `report_reproducibility_manifest.json`
- `report_reproducibility_checksums.tsv`
- `report_reproducibility_report.md`

## Scope

Checksums are computed over these generated reporting artifacts:
- `report_bundle_summary.md`
- `report_bundle_metrics.json`
- `report_bundle_top_signals.tsv`
- `report_dashboard.html`
- `report_dashboard_data.json`

## Checksum algorithm

- `sha256`

## Suggested usage

1. Run `--workflow report` on source environment.
2. Archive the three reproducibility artifacts with report outputs.
3. Re-run on target environment and compare `report_reproducibility_checksums.tsv`.
4. Investigate any hash drift before final release handoff.
