# Phase 5: RNA-seq Integration - Research

**Researched:** 2026-02-26
**Domain:** RNA-seq ingestion, dual normalization, and expression association modeling
**Confidence:** MEDIUM

<user_constraints>
## User Constraints (from 05-CONTEXT.md)

### Locked Decisions
- RNA source baseline: MAGE-linked public RNA-seq context.
- Emit both normalization outputs: log2(TPM+1) and quantile-normalized expression.
- Baseline association model: linear model on normalized expression.
- Baseline covariates: PC1-PC5, sex, batch, library metrics.

### Discretion
- Exact RNA metadata schema and mapping strategy.
- Library metrics field defaults and fallback policy.
- Output layout for dual normalization tracks.

### Deferred (OUT OF SCOPE)
- Mixed/random-effects RNA models.
- Final figure/dashboard packaging.
</user_constraints>

<research_summary>
## Summary

A robust Phase 5 baseline should separate ingestion/normalization from association execution while preserving sample-ID traceability across CN, variants, and RNA layers.

Dual normalization is best implemented as parallel artifacts built from a shared canonical expression matrix to avoid divergence and simplify downstream comparisons.

**Primary recommendation:** implement RNA ingestion + QC first, then run CN→expression and variant→expression models with identical covariate policy and standardized result schema.
</research_summary>

<standard_stack>
## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python (pandas/numpy) | current stable | RNA matrix transforms and joins | Efficient tabular handling |
| statsmodels | current stable | Linear model fitting | Interpretable inferential outputs |
| Nextflow DSL2 | current stable | Workflow orchestration | Portable and reproducible execution |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scipy | current stable | Quantile normalization helper ops | Needed for stable rank transforms |
| pyarrow/parquet | current stable | Efficient matrix serialization | For larger matrix scaling as needed |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Linear model baseline | GLM for all analyses | Flexible but higher complexity for baseline phase |
| Dual-output from one matrix | Separate independent pipelines | Easier to drift and harder to reconcile |
</standard_stack>

<architecture_patterns>
## Architecture Patterns

### Pattern 1: Canonical Matrix First
Create one canonical RNA expression matrix + metadata map before normalization branches.

### Pattern 2: Dual Normalization Branches
Compute `log2(TPM+1)` and quantile-normalized matrices from the canonical matrix with shared sample order and gene IDs.

### Pattern 3: Unified Association Schema
Emit CN-expression and variant-expression results with the same inferential field names where possible.

### Anti-Patterns to Avoid
- Inconsistent sample IDs across RNA and WGS-derived artifacts.
- Hidden normalization differences between CN and variant association paths.
- Missing technical covariate handling in linear models.
</architecture_patterns>

<common_pitfalls>
## Common Pitfalls

### Pitfall 1: Sample overlap collapse
Mitigate by explicit overlap reports before each model stage.

### Pitfall 2: Library metric missingness
Mitigate by defining required vs optional technical fields and fallback values with warnings.

### Pitfall 3: Normalization ambiguity
Mitigate by writing normalization metadata artifacts and clear output contracts.
</common_pitfalls>

<open_questions>
## Open Questions

1. **Gene set scope for EBV expression**
   - Recommendation: start with full available EBV gene features and allow optional subset lists.

2. **Quantile normalization tie handling**
   - Recommendation: deterministic rank averaging implementation with metadata snapshot.

3. **Covariate fallback behavior for missing library metrics**
   - Recommendation: explicit imputation strategy + warning counts in QC report.
</open_questions>

<sources>
## Sources

### Primary (HIGH confidence)
- Project planning artifacts and completed Phases 1-4 outputs.

### Secondary (MEDIUM confidence)
- Standard RNA-seq normalization and association workflow patterns.
</sources>
