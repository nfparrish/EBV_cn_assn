# Phase 4: WGS Association Modeling - Research

**Researched:** 2026-02-26
**Domain:** Covariate-adjusted variant association against EBV copy number
**Confidence:** MEDIUM

<user_constraints>
## User Constraints (from 04-CONTEXT.md)

### Locked Decisions
- Model family baseline: GLM.
- Outcome representation baseline: z-scored CN.
- Multiple testing baseline: FDR (Benjamini-Hochberg).
- Covariate policy: auto fallback for lab/batch handling.
- Maintain covariate intent: PC1-PC5, sex, and sequencing metadata.

### Discretion
- GLM family/link defaults and override mechanics.
- Chunking/execution strategy for scalable variant-wise modeling.
- Diagnostic output breadth.

### Deferred (OUT OF SCOPE)
- Mixed-model/kinship methods.
- RNA-seq modeling.
- Publication/dashboard packaging.
</user_constraints>

<research_summary>
## Summary

A practical Phase 4 baseline is a reproducible variant-wise GLM pipeline that standardizes CN outcome (z-score), harmonizes covariates, and emits per-variant statistics with FDR correction. This maintains interpretability while scaling to larger cohorts.

The safest implementation pattern is:
1) harmonize covariates and build model matrix,
2) run EBV and host variant association loops with consistent formulas,
3) emit diagnostics and adjusted results for QC/sign-off.

**Primary recommendation:** start with gaussian/identity GLM baseline for continuous z-scored CN and parameterize model family/link for later extension.
</research_summary>

<standard_stack>
## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python (statsmodels) | current stable | GLM fitting and model diagnostics | Mature and transparent statistical APIs |
| pandas/numpy | current stable | Matrix harmonization and table operations | Efficient cohort-scale data handling |
| Nextflow DSL2 | current stable | Orchestration and portability | Reproducible HPC execution |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| scipy | current stable | Statistical helper utilities | Additional tests/diagnostics |
| matplotlib/seaborn | current stable | Basic model QC plots | Diagnostic artifact generation |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| statsmodels GLM | sklearn + custom stats | Faster pipelines but weaker inferential/reporting transparency |
| Single-pass giant table | Chunked per-variant execution | Chunking improves memory behavior and parallelism |
</standard_stack>

<architecture_patterns>
## Architecture Patterns

### Pattern 1: Covariate Contract First
Build and validate a covariate matrix once, then join with each variant track.

### Pattern 2: Unified Model Interface
Use identical output schema for EBV and host variant association runs.

### Pattern 3: Built-In Multiple Testing
Run FDR adjustment as a required downstream stage, not an optional post-processing step.

### Anti-Patterns to Avoid
- Running models with silently missing covariates.
- Mixing different formulas across EBV/host tracks without provenance.
- Emitting p-values without adjusted p-values and diagnostics context.
</architecture_patterns>

<common_pitfalls>
## Common Pitfalls

### Pitfall 1: Sample mismatch across CN/variant/covariate tables
Mitigate with explicit overlap reporting before modeling.

### Pitfall 2: Outcome scale instability
Mitigate by deterministic z-score transformation with logged mean/sd.

### Pitfall 3: Over-interpretation without diagnostics
Mitigate by requiring model-level QC summaries and residual checks in outputs.
</common_pitfalls>

<open_questions>
## Open Questions

1. **Default GLM family/link**
   - Recommendation: gaussian/identity for v1, parameterized.

2. **Variant batching granularity**
   - Recommendation: chunk by variant blocks to balance overhead and memory.

3. **Minimum sample overlap threshold**
   - Recommendation: enforce a configurable minimum before fitting each variant model.
</open_questions>

<sources>
## Sources

### Primary (HIGH confidence)
- `.planning` requirements and completed Phase 2/3 outputs.

### Secondary (MEDIUM confidence)
- Standard large-scale GLM association workflow practices.
</sources>
