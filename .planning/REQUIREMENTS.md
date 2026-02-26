# Requirements: EBV Copy Number Multi-Omics Pipeline

**Defined:** 2026-02-26
**Core Value:** Deliver a reproducible, scheduler-agnostic pipeline that can process thousands of public samples and produce statistically robust EBV copy number association results.

## v1 Requirements

### Data Ingestion & QC

- [ ] **DATA-01**: Pipeline ingests 1000 Genomes Phase 3 WGS sample inputs from user-provided manifest or equivalent metadata table.
- [ ] **DATA-02**: Pipeline validates input integrity and writes standardized QC summary reports per sample and cohort.
- [ ] **DATA-03**: Pipeline supports resumable processing with explicit failed/complete sample status tracking.

### EBV Copy Number Estimation

- [ ] **CN-01**: Pipeline computes per-sample EBV copy number from WGS data using a documented and reproducible method.
- [ ] **CN-02**: Pipeline outputs normalized EBV copy number table with required metadata fields for downstream modeling.
- [ ] **CN-03**: Pipeline records confidence/quality metrics for copy-number estimates.

### Variant Calling

- [ ] **VAR-01**: Pipeline calls EBV genetic variants from WGS reads and exports cohort-level EBV variant matrix.
- [ ] **VAR-02**: Pipeline calls host (human) variants or consumes validated host variant input for association testing.
- [ ] **VAR-03**: Pipeline applies configurable filtering thresholds and logs retained/excluded variant counts.

### Association Analysis (WGS)

- [ ] **ASSOC-01**: Pipeline fits linear models for EBV copy number associations with covariates (PC1-PC5, sex, sequencing lab/batch).
- [ ] **ASSOC-02**: Pipeline outputs association statistics (effect size, standard error, p-value, adjusted p-value) for EBV and host variant analyses.
- [ ] **ASSOC-03**: Pipeline provides model diagnostics and summary plots for QC/interpretation.

### RNA-seq Expression Association

- [ ] **RNA-01**: Pipeline ingests selected public RNA-seq data source(s) for EBV expression analysis.
- [ ] **RNA-02**: Pipeline tests association between EBV copy number and EBV gene expression with documented covariates.
- [ ] **RNA-03**: Pipeline tests association between EBV genetic variants and EBV gene expression.

### Reporting & Portability

- [ ] **REP-01**: Pipeline generates manuscript-ready figures and publication-oriented result tables.
- [ ] **REP-02**: Pipeline generates an interactive dashboard for exploratory result review.
- [ ] **PORT-01**: Pipeline includes scheduler-agnostic execution profiles and setup documentation for HPC transfer.

## v2 Requirements

### Method Expansion

- **METH-01**: Add mixed-model support for population structure/relatedness beyond baseline linear models.
- **METH-02**: Add rare-variant burden analyses.

### Data Expansion

- **DATA-04**: Add additional large public WGS/RNA-seq cohorts beyond initial v1 data sources.

## Out of Scope

| Feature | Reason |
|---------|--------|
| Clinical-grade diagnostics and regulatory reporting | Not required for initial research pipeline goals |
| Automated private/protected dataset onboarding | v1 focuses on publicly available data |
| End-user no-code web productization | v1 prioritizes scientific pipeline correctness and portability |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| DATA-01 | Phase 1 | Pending |
| DATA-02 | Phase 1 | Pending |
| DATA-03 | Phase 1 | Pending |
| CN-01 | Phase 2 | Pending |
| CN-02 | Phase 2 | Pending |
| CN-03 | Phase 2 | Pending |
| VAR-01 | Phase 3 | Pending |
| VAR-02 | Phase 3 | Pending |
| VAR-03 | Phase 3 | Pending |
| ASSOC-01 | Phase 4 | Pending |
| ASSOC-02 | Phase 4 | Pending |
| ASSOC-03 | Phase 4 | Pending |
| RNA-01 | Phase 5 | Pending |
| RNA-02 | Phase 5 | Pending |
| RNA-03 | Phase 5 | Pending |
| REP-01 | Phase 6 | Pending |
| REP-02 | Phase 6 | Pending |
| PORT-01 | Phase 6 | Pending |

**Coverage:**
- v1 requirements: 18 total
- Mapped to phases: 18
- Unmapped: 0 ✓

---
*Requirements defined: 2026-02-26*
*Last updated: 2026-02-26 after initial definition*
