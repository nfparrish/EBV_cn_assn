---
phase: 01-data-intake-cohort-qc
plan: "02"
subsystem: workflow-orchestration
tags:
  - nextflow
  - ingest
  - resumability
provides:
  - Sample-level intake workflow
  - Per-sample status artifacts and run-level status registry
affects: []
tech-stack:
  added:
    - WGS_INTAKE DSL2 workflow
    - STAGE_INPUTS module
    - SAMPLE_STATUS module
  patterns:
    - Deterministic rerun behavior via status registry
key-files:
  created:
    - pipeline/main.nf
    - pipeline/nextflow.config
    - pipeline/workflows/wgs_intake.nf
    - pipeline/modules/local/stage_inputs/main.nf
    - pipeline/modules/local/sample_status/main.nf
  modified:
    - pipeline/conf/base.config
key-decisions:
  - decision: Status registry emitted as TSV plus JSON aggregate summary.
    rationale: TSV is easy to inspect and diff; JSON summary supports downstream tooling.
  - decision: Use optional status registry lookup to mark sample as skipped on rerun.
    rationale: Maintains deterministic resume while preserving visibility.
patterns-established:
  - Wave 1 execution supports parallel-ready sample handling.
duration: 65min
completed: 2026-02-26
---

# Phase 1: Data Intake & Cohort QC Summary

Implemented resumable sample intake workflow and status tracking foundations.

## Performance
- **Duration:** 65min
- **Tasks:** 2 completed
- **Files modified:** 6

## Accomplishments
- Added top-level Nextflow entrypoint and workflow wiring for intake execution.
- Implemented per-sample staging metadata generation with source existence checks.
- Implemented per-sample status output and cohort-level status registry/summary generation.

## Files Created/Modified
- `pipeline/main.nf` - Pipeline entrypoint.
- `pipeline/nextflow.config` - DSL2 config include.
- `pipeline/workflows/wgs_intake.nf` - Intake workflow and status summary process.
- `pipeline/modules/local/stage_inputs/main.nf` - Per-sample staging artifact generator.
- `pipeline/modules/local/sample_status/main.nf` - Resume-aware status emitter.
- `pipeline/conf/base.config` - Executor/resource profile defaults.

## Decisions & Deviations
- No functional deviation from plan objective.

## Next Phase Readiness
- Wave 2 (`01-03`) can consume sample status and staging outputs to produce QC report bundles.
