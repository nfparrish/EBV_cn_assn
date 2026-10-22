#!/usr/bin/env bash
#SBATCH --job-name=bcr_evid
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-731%20
#SBATCH --output=/work/nfp8/EBV/logs/evid_%A_%a.log
set -euo pipefail
BASE=/work/nfp8/EBV
sample=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1}' "$BASE/metadata/bcr_manifest.tsv")
[[ -n "${sample:-}" ]] || exit 1
out=$BASE/bcr_evidence/$sample.tsv
[[ -f $out ]] && exit 0
[[ -f $BASE/bcr/$sample/${sample}_report.tsv ]] || { echo "$sample: TRUST4 not done yet"; exit 0; }
mkdir -p "$BASE/bcr_evidence"
bash "$BASE/scripts/bcr_evidence.sh" "$sample" > "$out.tmp" && mv "$out.tmp" "$out"
cat "$out"
