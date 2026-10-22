#!/usr/bin/env bash
#SBATCH --job-name=mage_trust4
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-779%12
#SBATCH --output=/work/nfp8/EBV/logs/trust4_%A_%a.log
#
# BCR repertoire from MAGE LCL RNA-seq, via TRUST4 (assembles CDR3s from bulk RNA-seq,
# including partially-mapped and unmapped reads).
#
# Purpose: LCLs start polyclonal and become oligoclonal, so the repertoire is a per-line
# CULTURE-HISTORY probe -- the covariate IGSR and Coriell do not publish. Specifically:
#   - oligoclonality as a culture-age proxy;
#   - clone identity compared against the 2019 WGS and ~2023 ONT, giving a per-sample
#     measure of how far the cultures behind our predictor (CN) and outcome (lytic score)
#     have drifted apart -- which bounds H1 independent of biology;
#   - the 24 donors with 2 technical (same thaw) + 1 biological (different thaw) replicate
#     calibrate measurement noise against real thaw-to-thaw clonal turnover.
#
# Runs against whatever has finished downloading: a run without its .done marker is
# skipped, and resubmitting picks up newly completed runs. So this proceeds in parallel
# with the ongoing FASTQ download rather than waiting for it.
set -euo pipefail
BASE=/work/nfp8/MAGE
OUTD=/work/nfp8/EBV/trust4
T4=/hpc/group/parrishlab/envs/vcall-env/bin/run-trust4
REF=/hpc/group/parrishlab/refs/trust4

run=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$BASE/metadata/download_plan.tsv")
[[ -n "${run:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
[[ -f $BASE/fastq/$run/.done ]] || { echo "$run not downloaded yet -- skipping"; exit 0; }
[[ -f $OUTD/$run/${run}_report.tsv ]] && { echo "$run already processed"; exit 0; }

mkdir -p "$OUTD/$run"
r1=$BASE/fastq/$run/${run}_1.fastq.gz
r2=$BASE/fastq/$run/${run}_2.fastq.gz
echo "$(date -Is) $run on $(hostname)"

cd "$OUTD/$run"
$T4 -f "$REF/hg38_bcrtcr.fa" --ref "$REF/human_IMGT+C.fa" \
    -1 "$r1" -2 "$r2" -o "$run" -t "${SLURM_CPUS_PER_TASK:-4}" --stage 0

rep=$OUTD/$run/${run}_report.tsv
if [[ -f $rep ]]; then
  n=$(tail -n +2 "$rep" | wc -l)
  igh=$(awk -F'\t' 'NR>1 && $NF ~ /^IGH/ {c++} END{print c+0}' "$rep")
  top=$(awk -F'\t' 'NR>1 {s+=$1; if ($1>m) m=$1} END{printf "%.3f", (s>0 ? m/s : 0)}' "$rep")
  echo "  clonotypes=$n  IGH=$igh  top-clone fraction=$top"
  printf '%s\t%s\t%s\t%s\n' "$run" "$n" "$igh" "$top" >> /work/nfp8/EBV/trust4_summary.tsv
else
  echo "  WARNING: no report produced for $run"; exit 1
fi
