#!/usr/bin/env bash
# ==============================================================================
# Turn-Key HPC Slurm Orchestration Pipeline for Duke Compute Cluster (DCC)
# ==============================================================================
# Reproduces all raw-BAM/FASTQ alignments, variant calls, BCR assemblies, and
# consolidated master tables in `/work/nfp8/EBV/metadata/` (mirrored to
# `/hpc/group/parrishlab/EBV_1kGP_MAGE/data/` and `replication_2026/data/`).
#
# Usage:
#   # Mode 1: Fast recomputation of all statistical models & tables from
#   #         existing consolidated master tables (~30 seconds on login/compute node):
#   bash replication_2026/scripts/run_dcc_pipeline.sh --from-master-tables
#
#   # Mode 2: Full Slurm array submission from raw 1kGP CRAMs, ONT BAMs, and MAGE FASTQs:
#   bash replication_2026/scripts/run_dcc_pipeline.sh --submit-slurm
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EBV_BASE="${EBV_BASE:-/work/nfp8/EBV}"
PYTHON_BIN="${PYTHON_BIN:-/usr/bin/python3}"
PYSAM_PYTHON="${PYSAM_PYTHON:-/hpc/group/parrishlab/envs/pysam/bin/python3}"

MODE="${1:---from-master-tables}"

if [[ "${MODE}" == "--submit-slurm" ]]; then
    echo "[1/6] Submitting 1kGP 30x WGS, Phase 3 Low-Coverage, and ONT EBV CN + IR1 + GC-Bias arrays..."
    sbatch "${SCRIPT_DIR}/ebv_cn_slice.sbatch"
    sbatch "${SCRIPT_DIR}/ebv_cn_lowcov.sbatch"
    sbatch "${SCRIPT_DIR}/ir1_depth.sbatch"
    sbatch "${SCRIPT_DIR}/viral_gcbias.sbatch"
    sbatch "${SCRIPT_DIR}/ont_ir1.sbatch"

    echo "[2/6] Submitting LoFreq intrahost viral variant calling & oriLyt replication gradient check..."
    sbatch "${SCRIPT_DIR}/lofreq_call.sbatch"
    sbatch "${SCRIPT_DIR}/oriLyt_check.sbatch"

    echo "[3/6] Submitting WGS & MAGE RNA TRUST4 BCR repertoire assembly..."
    sbatch "${SCRIPT_DIR}/bcr_wgs.sbatch"
    sbatch "${SCRIPT_DIR}/trust4.sbatch"

    echo "[4/6] Submitting MAGE 779-library competitive splice-aware EBV RNA-seq quantification..."
    sbatch "${SCRIPT_DIR}/mage_ebv_rna.sbatch"

    echo "[5/6] Submitting within-sample heteroplasmic viral cis-ASE & fine-grained splice isoform array..."
    sbatch "${SCRIPT_DIR}/viral_ase.sbatch"

    echo "[6/6] Slurm jobs submitted. Once complete, run with --from-master-tables to regenerate all analysis tables."
    exit 0
fi

echo "=============================================================================="
echo "Running Turn-Key Analytical & Statistical Pipeline on Master Tables"
echo "=============================================================================="

if [[ -d "${EBV_BASE}/metadata" ]]; then
    echo "[Step 1] Running Trio Heritability Decomposition (heritability.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/heritability.py"

    echo "[Step 2] Running IR1 Population & GC-Controlled Analyses (ir1_pops.py, ir1_decisive.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/ir1_pops.py"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/ir1_decisive.py"

    echo "[Step 3] Running PoN-Cleaned BCR Clonality & RNA Recovery Analyses (clonality_vs_cn.py, rna_vs_wgs_clones.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/clonality_vs_cn.py"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/rna_vs_wgs_clones.py"

    echo "[Step 4] Running Full 779/779 MAGE H1 & Replicate Thaw Analysis (ebv_rna_analyse.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/ebv_rna_analyse.py"

    echo "[Step 5] Running Aim 1: All 84 Viral Genes & Abortive vs Productive Dissection (viral_gene_dissection.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/viral_gene_dissection.py"

    echo "[Step 6] Running Aim 2: Within-Population Viral CN/eQTLs & IR1 Structural QTL (viral_genetics_qtl.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/viral_genetics_qtl.py"

    echo "[Step 7] Consolidating Viral cis-ASE, Splice Isoforms, & Annotated LCL Catalog (consolidate_ase_and_ko.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/consolidate_ase_and_ko.py"

    echo "[Step 8] Running Aim 3: Baron-Kenny Host Transcriptome Mediation & Bipartite Map (host_mediation_bipartite.py)..."
    "${PYTHON_BIN}" "${SCRIPT_DIR}/host_mediation_bipartite.py"
fi

echo "[Step 9] Running Local Self-Contained Verification & Figure Generation (reproduce_all.py)..."
if command -v uv >/dev/null 2>&1; then
    uv run "${REPO_DIR}/reproduce_all.py"
else
    "${PYTHON_BIN}" "${REPO_DIR}/reproduce_all.py"
fi
