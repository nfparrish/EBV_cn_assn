# EBV_cn_assn

[![Docs](https://img.shields.io/badge/docs-index-blue)](pipeline/docs/INDEX.md)
[![HPC WGS10](https://img.shields.io/badge/HPC-WGS10%20sbatch-green)](pipeline/scripts/dcc_wgs10_testset.sbatch)

End-to-end Nextflow DSL2 pipeline for EBV copy-number and association analysis across WGS and RNA-seq, including reporting, dashboard generation, and reproducibility artifacts.

## HPC quick start (Duke DCC)

- Run 10-sample WGS test set via sbatch wrapper: [pipeline/scripts/dcc_wgs10_testset.sbatch](pipeline/scripts/dcc_wgs10_testset.sbatch)
- Usage and submit options: [pipeline/docs/hpc-portability.md](pipeline/docs/hpc-portability.md)

## Status

This repository contains the complete Phase 1–6 implementation:

- Phase 1: Data intake + cohort QC
- Phase 2: EBV copy number engine
- Phase 3: Variant calling + filtering/provenance
- Phase 4: WGS association modeling
- Phase 5: RNA integration (matrix and FASTQ recount modes)
- Phase 6: Reporting bundle + dashboard + reproducibility + HPC portability docs

Project closeout artifacts:
- `.planning/ROADMAP.md`
- `.planning/STATE.md`
- `.planning/COMPLETION-SUMMARY.md`

## Repository layout

- `pipeline/main.nf` — top-level workflow router
- `pipeline/workflows/` — workflow-level orchestration
- `pipeline/modules/local/` — DSL2 process modules
- `pipeline/scripts/` — Python processing/statistics/reporting scripts
- `pipeline/conf/base.config` — parameters, labels, profiles (`standard`, `slurm`, `pbs`)
- `pipeline/docs/` — contracts and operational documentation
- `pipeline/tests/fixtures/` — small fixtures/manifests for smoke runs

## Workflow entrypoints

Run from repo root:

- Intake:  
  `nextflow run pipeline/main.nf --workflow intake --manifest <manifest.csv> -profile standard`
- CN engine:  
  `nextflow run pipeline/main.nf --workflow cn --manifest <manifest.csv> -profile standard`
- Variant layer:  
  `nextflow run pipeline/main.nf --workflow variant --manifest <manifest.csv> --host_variant_input <host_variants.tsv> -profile standard`
- Association:  
  `nextflow run pipeline/main.nf --workflow assoc --manifest <manifest.csv> --host_variant_input <host_variants.tsv> -profile standard`
- RNA integration:  
  `nextflow run pipeline/main.nf --workflow rna --manifest <manifest.csv> -profile standard`
- Reporting bundle/dashboard/reproducibility:  
  `nextflow run pipeline/main.nf --workflow report --manifest <manifest.csv> --host_variant_input <host_variants.tsv> -profile standard`

## RNA modes

Configured via `--rna_input_mode`:

- `matrix` — consume expression matrix input/fallback generation
- `fastq` — align/count from FASTQ against combined host+EBV reference
- `auto` — choose FASTQ mode when references are supplied and matrix input is absent

FASTQ mode parameters:

- `--rna_reference_fasta <combined_host_ebv.fa>`
- `--rna_reference_gtf <combined_host_ebv.gtf>`
- optional: `--rna_aligner minimap2`, `--rna_counter auto`, `--rna_threads <n>`

## Reporting outputs

`--workflow report` emits:

- `report_bundle_summary.md`
- `report_bundle_metrics.json`
- `report_bundle_top_signals.tsv`
- `report_dashboard.html`
- `report_dashboard_data.json`
- `report_reproducibility_manifest.json`
- `report_reproducibility_checksums.tsv`
- `report_reproducibility_report.md`

## Core software inventory

Tracked core implementation files currently include:

- 6 workflow files in `pipeline/workflows/`
- 33 module files in `pipeline/modules/local/*/main.nf`
- 29 Python scripts in `pipeline/scripts/`

## Documentation index

- Canonical docs index: `pipeline/docs/INDEX.md`

- Manifest spec: `pipeline/docs/manifest-spec.md`
- CN outputs: `pipeline/docs/ebv-cn-output-contract.md`
- Variant layer: `pipeline/docs/ebv-variant-calling.md`, `pipeline/docs/variant-filtering.md`
- Association outputs: `pipeline/docs/association-output-contract.md`, `pipeline/docs/association-diagnostics.md`
- RNA contracts: `pipeline/docs/rna-input-normalization.md`, `pipeline/docs/rna-association-contract.md`, `pipeline/docs/rna-association-qc.md`, `pipeline/docs/rna-variant-presence-qc.md`, `pipeline/docs/rna-fastq-recount.md`
- Reporting/dashboard: `pipeline/docs/reporting-bundle.md`, `pipeline/docs/report-dashboard.md`
- Portability/reproducibility: `pipeline/docs/hpc-portability.md`, `pipeline/docs/reproducibility-checks.md`

## Notes

- Large local test datasets are intentionally excluded from git (`pipeline/tests/data/` is ignored).
- Nextflow `work/` and local runtime logs are excluded from version control.
