# EBV CN confidence scoring

## Overview

Phase 2 Wave 2 adds rule-based confidence annotations for each sample-level EBV CN estimate and cohort-level QC outputs to support pilot sign-off.

Inputs:

- `ebv_cn_normalized.tsv`
- `ebv_cn_cohort_summary.json`

Outputs:

- `ebv_cn_confidence.tsv`
- `ebv_cn_confidence_summary.json`
- `ebv_cn_qc_summary.json`
- `ebv_cn_qc_report.md`
- `ebv_cn_low_confidence.tsv`

## Rule-based confidence model

A per-sample score starts at `100` and subtracts penalties:

- `-20` if `depth_ratio_cn < 0.5` (`low_depth_signal`)
- `-20` if `read_count_ratio_cn < 1.0` (`low_readcount_signal`)
- `-30` if `abs(depth_ratio_cn_norm - read_count_ratio_cn_norm) > 0.5` (`method_discordance_high`)
- `-15` if `abs(depth_ratio_cn_norm - read_count_ratio_cn_norm) > 0.25` (`method_discordance_moderate`)
- `-20` if either normalized method value `< 0.25` (`extreme_normalized_signal`)

Score is clipped to `[0, 100]`.

### Confidence tiers

- `high`: score `>= 80`
- `medium`: score `>= 55` and `< 80`
- `low`: score `< 55`

## QC interpretation

`ebv_cn_qc_summary.json` reports:

- score mean/median
- tier counts
- low-confidence count and rate
- per-flag counts

`ebv_cn_low_confidence.tsv` is the review list for manual inspection and threshold tuning.

## Notes

- This Wave 2 implementation is explicitly rule-based and deterministic.
- The design is intended to be replaced or extended by model-based confidence calibration in later phases if required.
