#!/usr/bin/env bash
#SBATCH --job-name=ebv_rna_scav
#SBATCH --account=parrishlab
#SBATCH --partition=scavenger
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --array=1-669%120
#SBATCH --output=/work/nfp8/EBV/logs/ebv_rna_scav_%A_%a.log
set -euo pipefail

orig_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /work/nfp8/EBV/metadata/scav_order.txt)
[[ -n "${orig_id:-}" ]] || exit 0
SLURM_ARRAY_TASK_ID=$orig_id exec /work/nfp8/EBV/scripts/mage_ebv_rna.sbatch
