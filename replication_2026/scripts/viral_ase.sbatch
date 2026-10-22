#!/usr/bin/env bash
#SBATCH --job-name=viral_ase
#SBATCH --account=parrishlab
#SBATCH --partition=common,scavenger
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --array=1-60%60
#SBATCH --output=/work/nfp8/EBV/logs/vase_%A_%a.log
set -euo pipefail

BASE=/work/nfp8/MAGE
OUTD=/work/nfp8/EBV/ase_mage
MM2=/hpc/group/parrishlab/envs/vcall-env/bin/minimap2
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
EBV_REF=/hpc/group/parrishlab/refs/ebv/chrEBV.fa

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /work/nfp8/EBV/metadata/ase_60_runs.txt)
[[ -n "${run:-}" ]] || exit 0
[[ -s "$OUTD/${run}.ase.tsv" ]] && exit 0

r1=$BASE/fastq/$run/${run}_1.fastq.gz
r2=$BASE/fastq/$run/${run}_2.fastq.gz
tmp=$(mktemp -d "${TMPDIR:-/tmp}/vase_${run}_XXXXXX")
trap 'rm -rf "$tmp"' EXIT

$MM2 -t 2 -ax splice:sr --secondary=no "$EBV_REF" "$r1" "$r2" 2>/dev/null \
  | $SAM view -h -F 2308 - > "$tmp/ebv.sam"

python3 /work/nfp8/EBV/scripts/viral_ase_worker.py "$run" "$tmp/ebv.sam" "$OUTD/${run}.ase.tsv" "$OUTD/${run}.spl.tsv"
