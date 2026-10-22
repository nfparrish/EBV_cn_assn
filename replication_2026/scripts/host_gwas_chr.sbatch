#!/usr/bin/env bash
#SBATCH --job-name=ebv_hgwas
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --array=1-22
#SBATCH --output=/work/nfp8/EBV/logs/hgwas_%A_%a.log
set -euo pipefail

PY311=/hpc/group/parrishlab/ms1299/miniconda3/envs/py311_bash5_2/bin/python3
"$PY311" /work/nfp8/EBV/scripts/run_host_gwas_chr.py "$SLURM_ARRAY_TASK_ID"
