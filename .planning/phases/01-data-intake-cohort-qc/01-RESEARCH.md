# Phase 1: Data Intake & Cohort QC - Research

**Researched:** 2026-02-26
**Domain:** Large-scale WGS ingest and QC for HPC-portable workflow systems
**Confidence:** MEDIUM

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Use Nextflow DSL2 with nf-core style module boundaries.
- 1000 Genomes Phase 3 is the required initial WGS source.
- Scheduler-agnostic HPC portability is required.
- Output per-sample and cohort-level QC artifacts for downstream phases.
- Support resumable processing that avoids rerunning successful work.

### Claude's Discretion
- Manifest schema specifics and optional metadata fields.
- Exact QC metric mix and threshold defaults.
- Concrete project directory layout as long as nf-core-compatible.
- Container strategy across local/HPC contexts.

### Deferred Ideas (OUT OF SCOPE)
- Association analysis implementation (Phase 4).
- RNA-seq expression workflow (Phase 5).
- Final dashboard/report polish (Phase 6).
- Mixed-model/rare-variant methods (v2).
</user_constraints>

<research_summary>
## Summary

Phase 1 should implement a robust intake-and-validate backbone before any EBV-specific analytics. For large public WGS cohorts, the highest-value pattern is to define a strict sample manifest contract, normalize metadata early, and persist machine-readable run manifests (success/fail/skip) that make resume behavior deterministic.

For HPC portability, the standard path is Nextflow DSL2 with profile-driven configuration, containerized processes, and explicit resource labels for scheduler translation. Using staged channels (manifest -> validate -> fetch/locate inputs -> QC metrics -> consolidated report) keeps data lineage clear and reduces phase-to-phase coupling.

**Primary recommendation:** Build Phase 1 as a manifest-first, profile-driven Nextflow workflow with strict validation gates and MultiQC-backed cohort reporting.
</research_summary>

<standard_stack>
## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Nextflow DSL2 | current stable | Workflow orchestration | Native HPC/cloud portability and resume semantics |
| nf-core layout patterns | current conventions | Pipeline structure and module boundaries | Proven conventions for reusable genomics workflows |
| MultiQC | current stable | Cohort-level QC aggregation | Widely adopted consolidated QC reporting |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| FastQC | current stable | Read-level QC snapshots | Input FASTQ-based ingest path |
| samtools stats/flagstat | current stable | Alignment-level QC metrics | BAM/CRAM quality summaries |
| mosdepth or depth utility | current stable | Coverage-oriented QC metrics | Depth/distribution checks for downstream CN readiness |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Nextflow | Snakemake | Python-native ergonomics but weaker built-in profile portability patterns for this team preference |
| MultiQC aggregation | Custom HTML scripts | More flexibility but unnecessary maintenance burden |

**Installation (implementation phase):**
```bash
# example baseline (exact tools decided during implementation)
nextflow run main.nf -profile test
```
</standard_stack>

<architecture_patterns>
## Architecture Patterns

### Recommended Project Structure
```
pipeline/
├── main.nf
├── nextflow.config
├── conf/
│   ├── base.config
│   ├── slurm.config
│   └── pbs.config
├── modules/
│   ├── local/
│   └── nf-core/
├── assets/
│   └── schema/
└── workflows/
```

### Pattern 1: Manifest-First Validation Gate
**What:** Validate schema and file accessibility before expensive QC starts.
**When to use:** Always, especially for public cohort ingest with heterogeneous metadata.

### Pattern 2: Deterministic Resume via Sample Status Registry
**What:** Emit per-sample state (`validated`, `processed`, `failed`) and use it on reruns.
**When to use:** Any cohort larger than pilot scale.

### Anti-Patterns to Avoid
- **Implicit metadata assumptions:** downstream failures become hard to debug; validate early with explicit reasons.
- **Cluster-specific hard-coding:** breaks portability and violates phase constraints.
- **One giant process block:** reduces observability and reuse; prefer modular DSL2 processes/workflows.
</architecture_patterns>

<dont_hand_roll>
## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Cohort QC dashboarding | Custom report assembler | MultiQC | Existing ecosystem solves data aggregation and visualization well |
| HPC resume bookkeeping | Ad hoc checkpoint files | Nextflow native caching + explicit status outputs | Avoid brittle state handling |
| Metadata validation parser from scratch | Loose grep/awk checks only | Structured schema validation (JSON schema + scripted checks) | Clear errors and maintainable contracts |

**Key insight:** Most Phase-1 risk is operational robustness, not algorithmic novelty; prioritize proven orchestration and reporting tools.
</dont_hand_roll>

<common_pitfalls>
## Common Pitfalls

### Pitfall 1: Manifest Drift Across Batches
**What goes wrong:** Different batches include inconsistent columns or naming conventions.
**How to avoid:** Versioned schema + strict validation + normalization layer.

### Pitfall 2: Silent Partial Failures
**What goes wrong:** Long HPC runs complete with hidden per-sample failures.
**How to avoid:** Mandatory per-sample status table and fail-summary artifact in every run.

### Pitfall 3: Non-portable Resource Profiles
**What goes wrong:** Pipeline tuned for one cluster and fails on another scheduler.
**How to avoid:** Base profile + scheduler-specific profile overrides with documented defaults.
</common_pitfalls>

<open_questions>
## Open Questions

1. **Exact pilot subset strategy for phase validation**
   - Recommendation: define a representative small manifest (coverage/diversity) for smoke tests.

2. **Preferred container backend per HPC environment**
   - Recommendation: support both Docker (local) and Singularity/Apptainer (HPC) profiles.

3. **Minimum QC thresholds for early exclusion**
   - Recommendation: start with descriptive metrics + warning flags; enforce hard filters in later iteration after pilot review.
</open_questions>

<sources>
## Sources

### Primary (HIGH confidence)
- Project constraints and roadmap documents in `.planning/`.
- Established Nextflow/nf-core/MultiQC usage conventions (to be pinned with exact version docs during implementation).

### Secondary (MEDIUM confidence)
- Common large-cohort WGS workflow patterns from public bioinformatics practice.

### Tertiary (LOW confidence)
- None.
</sources>
