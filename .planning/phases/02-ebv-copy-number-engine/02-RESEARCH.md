# Phase 2: EBV Copy Number Engine - Research

**Researched:** 2026-02-26
**Domain:** EBV copy number estimation from WGS with cohort normalization and confidence scoring
**Confidence:** MEDIUM

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Consume Phase 1 validated intake/QC outputs as upstream input.
- Compute per-sample EBV copy number.
- Emit normalized cohort-ready CN table with required metadata.
- Emit confidence and quality annotations per sample.
- Keep implementation portable via Nextflow DSL2 and profile-driven HPC execution.

### Claude's Discretion
- EBV reference indexing/packaging details.
- Exact depth-based CN formula details.
- Confidence-scoring formula and threshold defaults.
- Internal pipeline file organization.

### Deferred Ideas (OUT OF SCOPE)
- Association testing workflows (Phase 4+).
- RNA-seq analysis (Phase 5).
- Final reporting/dashboard refinement (Phase 6).
</user_constraints>

<research_summary>
## Summary

The safest Phase 2 implementation pattern is depth-based EBV CN estimation with explicit host-normalized scaling and robust metadata propagation. For WGS cohorts, reliable CN is more about calibration and quality controls than raw counting alone.

A practical architecture is: reference prep -> EBV-aligned depth extraction -> host baseline extraction -> CN estimation -> cohort normalization table -> confidence scoring artifacts. This keeps method assumptions visible and supports reproducibility across HPC targets.

**Primary recommendation:** Build a transparent depth-ratio CN workflow with explicit confidence flags and pilot calibration against known-quality samples.
</research_summary>

<standard_stack>
## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| samtools | current stable | BAM/CRAM region/stat extraction | Reliable low-level alignment metrics |
| bedtools/mosdepth (or equivalent) | current stable | Depth estimation across target regions | Common depth quantification pattern |
| Nextflow DSL2 | current stable | Portable orchestration | Reproducible HPC execution |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Python (pandas/numpy) | current stable | CN table assembly and normalization | Post-depth metric synthesis |
| Reference index tooling | current stable | FASTA/aux index preparation | Deterministic reference handling |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Depth-ratio CN baseline | k-mer based quick estimators | Faster in some contexts but less transparent for this workflow stage |

</standard_stack>

<architecture_patterns>
## Architecture Patterns

### Pattern 1: Reference-First Determinism
Prepare and version EBV reference assets once, then reuse across samples.

### Pattern 2: Dual-Signal CN Estimation
Use EBV signal + host baseline context to reduce spurious CN inflation from uneven coverage.

### Pattern 3: Confidence Tiering
Emit confidence labels (`high`, `medium`, `low`) based on depth support and QC warnings.

### Anti-Patterns to Avoid
- Hidden normalization assumptions not encoded in outputs.
- CN outputs without sample-level quality/context fields.
- Coupling CN computation directly to downstream association code.
</architecture_patterns>

<common_pitfalls>
## Common Pitfalls

### Pitfall 1: Coverage Heterogeneity Misread as CN Signal
Mitigate by carrying Phase 1 QC metrics and explicit normalization factors.

### Pitfall 2: Incomplete Metadata Carry-Through
Mitigate by enforcing output schema that includes `sample_id`, `population`, `sex`, `sequence_lab`, QC flags.

### Pitfall 3: Uninterpretable Confidence Scores
Mitigate by documenting confidence formula and exposing constituent metrics.
</common_pitfalls>

<open_questions>
## Open Questions

1. **Reference target definition granularity**
   - Recommendation: start with a documented canonical EBV reference and fixed target windows.

2. **Normalization baseline strategy**
   - Recommendation: implement default host-depth scaling with override hooks.

3. **Confidence thresholds**
   - Recommendation: start with conservative defaults and tune after pilot review.
</open_questions>

<sources>
## Sources

### Primary (HIGH confidence)
- `.planning` requirements, roadmap, and Phase 1 outputs.

### Secondary (MEDIUM confidence)
- Common WGS depth-based CN practice patterns for viral copy number workflows.
</sources>
