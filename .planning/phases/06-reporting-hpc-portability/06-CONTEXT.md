---
phase: "06"
name: "Reporting & HPC Portability"
created: 2026-02-26
---

# Phase 6: Reporting & HPC Portability — Context

## Decisions

- Phase 6 consumes completed outputs from:
  - Phase 4 (`assoc` diagnostics/QC/top-hits)
  - Phase 5 (`rna` association QC/top-hits/variant presence)
- Wave 1 focuses on **reporting bundle generation** (figures/tables/narrative summary artifacts).
- Keep baseline implementation scheduler-agnostic and workflow-native in Nextflow DSL2.

## Discretion Areas

- Exact markdown narrative structure and summary metric schema.
- Top-hit table merge strategy across WGS and RNA layers.
- Naming of report bundle artifacts, as long as contracts are documented.

## Deferred to later waves

- Interactive dashboard packaging (Wave 2).
- Full HPC deployment packaging and reproducibility checklists (Wave 3).

## Phase Boundary

Phase 6 Wave 1 delivers automated report bundle artifacts suitable for manuscript-oriented review and handoff.
