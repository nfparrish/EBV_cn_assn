# Phase 3: Variant Calling Layer - Research

**Researched:** 2026-02-26
**Domain:** EBV + host variant generation with reproducible filtering and provenance
**Confidence:** MEDIUM

<user_constraints>
## User Constraints (from 03-CONTEXT.md)

### Locked Decisions
- EBV caller baseline is GATK Mutect2.
- Host variant path is consume public validated host VCF/PLINK.
- Filtering posture is balanced.
- Emit both AF and binary EBV matrices.
- Keep implementation portable in Nextflow DSL2.

### Discretion
- Balanced threshold exact defaults and tuning knobs.
- EBV matrix representation details beyond required outputs.
- Host harmonization internals for association-ready format.

### Deferred (OUT OF SCOPE)
- Association model execution.
- RNA-seq workflows.
- Final reporting/dashboard productization.
</user_constraints>

<research_summary>
## Summary

For Phase 3, the most robust architecture is a three-stage layer:
1) EBV variant calling + cohort aggregation,
2) host variant intake/harmonization,
3) balanced filtering with machine-readable provenance.

Mutect2 provides a strong baseline for EBV variant detection in mixed-depth contexts, while consuming public host variants reduces compute and avoids introducing an unnecessary host-calling dependency in v1.

**Primary recommendation:** implement a deterministic EBV call-to-matrix path, enforce cohort sample-ID harmonization against host variant inputs, and emit filtering provenance artifacts as first-class outputs.
</research_summary>

<standard_stack>
## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| GATK Mutect2 | current stable | EBV variant calling baseline | Mature caller with strong filtering ecosystem |
| bcftools | current stable | VCF normalization/index/query | Standard utility layer for matrix/provenance prep |
| Nextflow DSL2 | current stable | Portable orchestration | Reproducible HPC execution |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Python (pandas/pyarrow optional) | current stable | Variant matrix synthesis + provenance reports | Transform/aggregate VCF-derived tables |
| PLINK 2 | current stable | Host variant format harmonization | Standard host association input path |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Mutect2 baseline | bcftools-only calling | Simpler but may be less robust for low-fraction EBV variants |
| Public host-variant ingestion | In-pipeline host calling | More complete provenance but high compute and complexity |
</standard_stack>

<architecture_patterns>
## Architecture Patterns

### Pattern 1: Call → Normalize → Aggregate
Isolate per-sample EBV call artifacts before cohort matrix assembly to improve resumability.

### Pattern 2: Dual Matrix Contract
Emit both AF and binary matrices from the same canonical intermediate table to prevent drift.

### Pattern 3: Filter Ledger
Record retained/excluded counts by rule and stage in JSON/TSV artifacts for reproducibility.

### Anti-Patterns to Avoid
- Hard-coded thresholds without config logging.
- Mismatched sample IDs between EBV and host pathways.
- Only emitting final matrix without filter lineage.
</architecture_patterns>

<common_pitfalls>
## Common Pitfalls

### Pitfall 1: Variant normalization drift
Mitigate with explicit VCF normalization/order/indexing steps before matrix generation.

### Pitfall 2: Host input mismatch (build/sample IDs)
Mitigate with mandatory harmonization report and fail-fast mismatch counts.

### Pitfall 3: Non-reproducible filtering
Mitigate by writing threshold config snapshot and retained/excluded tables per stage.
</common_pitfalls>

<open_questions>
## Open Questions

1. **Mutect2 panel-of-normals usage in v1**
   - Recommendation: optional hook, off by default for pilot.

2. **Balanced default thresholds**
   - Recommendation: start with moderate QUAL/DP/AF defaults and expose params for calibration.

3. **Binary matrix threshold from AF**
   - Recommendation: configurable default cutoff with explicit reporting.
</open_questions>

<sources>
## Sources

### Primary (HIGH confidence)
- `.planning` requirements, roadmap, and completed Phase 1-2 outputs.

### Secondary (MEDIUM confidence)
- Common production patterns for Mutect2-based variant workflows and VCF harmonization.
</sources>
