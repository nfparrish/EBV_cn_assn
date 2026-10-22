# EBV_cn_assn

[![Docs](https://img.shields.io/badge/docs-index-blue)](pipeline/docs/INDEX.md)
[![HPC WGS10](https://img.shields.io/badge/HPC-WGS10%20sbatch-green)](pipeline/scripts/dcc_wgs10_testset.sbatch)
[![2026 Replication](https://img.shields.io/badge/Replication%202026-14%2F14%20Verified-brightgreen)](replication_2026/README.md)

End-to-end Nextflow DSL2 pipeline and 2026 full-cohort multi-omics analysis (`3,200` NYGC 30x WGS, `2,499` Phase 3 low-coverage WGS, `907` HPRC/1kGP ONT long-read WGS, and `779` MAGE RNA-seq LCLs) for Epstein-Barr virus (EBV) copy-number architecture, viral genetics, and host transcriptome mediation.

## 2026 Full-Cohort Replication & Multi-Omics Package (`replication_2026/`)

A self-contained turn-key replication bundle for all 2026 paper findings, master data tables, publication figures, and Duke Compute Cluster (DCC) Slurm orchestration is located in [`replication_2026/`](replication_2026/README.md):

```bash
# Verify all 14/14 headline statistical claims and regenerate all 6 publication figures (<15s)
uv run replication_2026/reproduce_all.py
```

- **Turn-Key Verifier**: [`replication_2026/reproduce_all.py`](replication_2026/reproduce_all.py) -> generates [`replication_2026/REPLICATION_VERIFICATION_REPORT.md`](replication_2026/REPLICATION_VERIFICATION_REPORT.md)
- **Figure Generator**: [`replication_2026/scripts/generate_figures.py`](replication_2026/scripts/generate_figures.py) -> outputs 6 multi-panel 300 DPI `.png` and `.pdf` figures in [`replication_2026/figures/`](replication_2026/figures/)
- **DCC HPC Orchestrator**: [`replication_2026/scripts/run_dcc_pipeline.sh`](replication_2026/scripts/run_dcc_pipeline.sh) (`--stage wgs|rnaseq|downstream|verify|all`)
- **Documentation & Claims Ledger**: [`replication_2026/docs/HANDOFF.md`](replication_2026/docs/HANDOFF.md) and [`replication_2026/docs/CLAIMS_LEDGER.md`](replication_2026/docs/CLAIMS_LEDGER.md)
- **Persistent DCC Mirror**: `/hpc/group/parrishlab/EBV_1kGP_MAGE/` (`nfp8@dcc-login.oit.duke.edu`)

---

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
- Heritability (parent-offspring trio analysis):  
  `nextflow run pipeline/main.nf --workflow heritability --manifest <pedigree_manifest.csv> --cram_reference_fasta <GRCh38.fa> -profile standard`

## Heritability analysis

The `heritability` workflow estimates h² for EBV copy number using parent-offspring trios.
Requires a manifest with pedigree columns (`family_id`, `father_id`, `mother_id`, `relationship`).

Methods (following Houldcroft et al. PLOS ONE 2014; Mandage et al. PLOS ONE 2017):
- Mid-parent regression (Falconer's method): h² = slope(offspring_CN ~ mid_parent_CN)
- Intraclass correlation (ICC) among family members
- Permutation test (10,000 shuffles) for empirical p-value
- Bootstrap (1,000 resamples) for 95% CI

Pedigree manifests can be generated from the 1000G `.ped` file:
```bash
python pipeline/scripts/build_pedigree_manifest.py manifests/1000g_pedigree.ped output.csv --max-trios 34
```

Key parameters:
- `--heritability_n_permutations` (default: 10000)
- `--heritability_n_bootstrap` (default: 1000)
- `--cram_reference_fasta` — required when using remote CRAMs (e.g., GRCh38 reference for 1000G 30x data)

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

## Trio manifests

- `manifests/trio_100_30x_manifest.csv` — 102 samples (34 trios) with 1000G 30x CRAM FTP URLs
- `manifests/all_trios_30x_manifest.csv` — 1910 samples (641 trios) for full heritability analysis
- `manifests/1000g_pedigree.ped` — 1000 Genomes pedigree file (3501 individuals)

## Notes

- Large local test datasets are intentionally excluded from git (`pipeline/tests/data/` is ignored).
- Nextflow `work/` and local runtime logs are excluded from version control.
- The `--cram_reference_fasta` parameter is required when processing remote CRAMs (samtools needs the reference for CRAM decoding).
