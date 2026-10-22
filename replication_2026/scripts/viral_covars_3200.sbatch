#!/usr/bin/env bash
#SBATCH --job-name=ebv_vcov3200
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --array=1-16
#SBATCH --output=/work/nfp8/EBV/logs/vcov3200_%A_%a.log
set -euo pipefail

PYSAM_PY=/hpc/group/parrishlab/envs/pysam/bin/python3
"$PYSAM_PY" /work/nfp8/EBV/scripts/viral_covars_3200.py "$SLURM_ARRAY_TASK_ID" 16
