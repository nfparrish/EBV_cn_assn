# Roadmap: EBV Copy Number Multi-Omics Pipeline

## Overview

This roadmap progresses from reliable data intake and EBV copy-number estimation to variant association analyses, RNA-seq integration, and publication-ready outputs, while ensuring the full workflow is portable to multiple HPC schedulers.

## Phases

- [x] **Phase 1: Data Intake & Cohort QC** - Build robust WGS ingest, manifest handling, and cohort/sample QC foundations.
- [x] **Phase 2: EBV Copy Number Engine** - Implement and validate per-sample EBV copy-number estimation.
- [x] **Phase 3: Variant Calling Layer** - Produce EBV and host variant datasets with filtering and provenance.
- [x] **Phase 4: WGS Association Modeling** - Run covariate-adjusted association analysis for EBV CN with variant features.
- [x] **Phase 5: RNA-seq Integration** - Connect EBV CN and EBV variants to EBV gene expression in public RNA-seq cohorts.
- [x] **Phase 6: Reporting & HPC Portability** - Deliver dashboard/figures/reports and scheduler-agnostic operational packaging.

## Phase Details

### Phase 1: Data Intake & Cohort QC
**Goal**: Reliable and resumable ingestion of public WGS inputs with baseline quality reporting.
**Depends on**: Nothing (first phase)
**Requirements**: DATA-01, DATA-02, DATA-03
**Success Criteria** (what must be TRUE):
  1. User can provide WGS sample manifest and run ingestion end-to-end.
  2. QC report artifacts are generated for each sample and summarized at cohort level.
  3. Failed samples are tracked and reruns only process incomplete work.
**Plans**: 3 plans

Plans:
- [x] 01-01: Define canonical manifests, metadata schema, and validation checks.
- [x] 01-02: Implement Nextflow ingestion and resumable sample-level execution.
- [x] 01-03: Add QC metrics computation and report generation.

### Phase 2: EBV Copy Number Engine
**Goal**: Produce trusted per-sample EBV copy-number values with quality/confidence annotations.
**Depends on**: Phase 1
**Requirements**: CN-01, CN-02, CN-03
**Success Criteria** (what must be TRUE):
  1. EBV copy number is computed for each processed sample.
  2. Output includes normalized cohort-level CN table suitable for modeling.
  3. Confidence/quality metrics accompany each CN estimate.
**Plans**: 3 plans

Plans:
- [x] 02-01: Implement EBV reference handling and copy-number algorithm process.
- [x] 02-02: Add normalization and metadata harmonization for CN outputs.
- [x] 02-03: Validate CN output quality on a pilot subset.

### Phase 3: Variant Calling Layer
**Goal**: Generate analysis-ready EBV and host variant datasets with filtering provenance.
**Depends on**: Phase 2
**Requirements**: VAR-01, VAR-02, VAR-03
**Success Criteria** (what must be TRUE):
  1. EBV variant matrix is generated for cohort samples.
  2. Host variant input is available in compatible format for association tests.
  3. Filtering thresholds and retained/excluded counts are reproducibly reported.
**Plans**: 3 plans

Plans:
- [x] 03-01: Implement EBV variant calling and aggregation outputs.
- [x] 03-02: Integrate host variant call/input pathway.
- [x] 03-03: Add filtering configuration and variant provenance reports.

### Phase 4: WGS Association Modeling
**Goal**: Deliver baseline linear association analyses linking variants and EBV copy number.
**Depends on**: Phase 3
**Requirements**: ASSOC-01, ASSOC-02, ASSOC-03
**Success Criteria** (what must be TRUE):
  1. Linear models run with required covariates (PC1-PC5, sex, sequencing lab/batch).
  2. Output tables include effect size, standard error, p-value, and adjusted p-value.
  3. Diagnostic summaries/plots are produced for model quality review.
**Plans**: 3 plans

Plans:
- [x] 04-01: Build covariate harmonization and model-ready matrices.
- [x] 04-02: Implement EBV and host variant association workflows.
- [x] 04-03: Add diagnostics and statistical QC artifacts.

### Phase 5: RNA-seq Integration
**Goal**: Quantify relationships between EBV CN/EBV variants and EBV gene expression.
**Depends on**: Phase 4
**Requirements**: RNA-01, RNA-02, RNA-03
**Success Criteria** (what must be TRUE):
  1. RNA-seq cohort(s) are ingested into a harmonized analysis format.
  2. EBV CN-to-expression associations are computed and reported.
  3. EBV variant-to-expression associations are computed and reported.
**Plans**: 3 plans

Plans:
- [x] 05-01: Define RNA-seq cohort ingestion and normalization pathway.
- [x] 05-02: Implement CN-expression association models.
- [x] 05-03: Implement EBV variant-expression association models.

### Phase 6: Reporting & HPC Portability
**Goal**: Produce publication-oriented outputs and make execution portable across HPC systems.
**Depends on**: Phase 5
**Requirements**: REP-01, REP-02, PORT-01
**Success Criteria** (what must be TRUE):
  1. Manuscript-ready figures and summary tables are generated from pipeline outputs.
  2. Interactive dashboard renders key outcomes for exploratory review.
  3. Workflow runs with scheduler profiles and documented setup on multiple HPC environments.
**Plans**: 3 plans

Plans:
- [x] 06-01: Implement reporting bundle (figures, tables, narrative summary).
- [x] 06-02: Build interactive dashboard output package.
- [x] 06-03: Finalize scheduler profiles, deployment docs, and reproducibility checks.

## Progress

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Data Intake & Cohort QC | 3/3 | Complete | 2026-02-26 |
| 2. EBV Copy Number Engine | 3/3 | Complete | 2026-02-26 |
| 3. Variant Calling Layer | 3/3 | Complete | 2026-02-26 |
| 4. WGS Association Modeling | 3/3 | Complete | 2026-02-26 |
| 5. RNA-seq Integration | 3/3 | Complete | 2026-02-26 |
| 6. Reporting & HPC Portability | 3/3 | Complete | 2026-02-26 |
