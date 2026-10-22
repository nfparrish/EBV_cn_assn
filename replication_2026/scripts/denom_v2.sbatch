#!/usr/bin/env bash
#SBATCH --job-name=ebv_denom
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=1-208%15
#SBATCH --output=/work/nfp8/EBV/logs/denom_%A_%a.log
#
# Robust autosomal denominator: per-window mean depth over 150 x 100 kb windows
# (15 Mb, evenly spread across the autosomes), rather than one pooled count over
# 24 windows. Saves every window so bias can be audited without re-fetching.
#
# `samtools depth` needs POS and CIGAR but not SEQ, so no human reference is needed
# even though this reads autosomal data.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
unset REF_PATH REF_CACHE
WINDOWS=$BASE/metadata/windows_v2.bed
MANIFEST=${MANIFEST:-$BASE/metadata/validation_subset.tsv}

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$MANIFEST")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/denom/$sample.windows.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/denom"
echo "$(date -Is) $sample on $(hostname)"

for a in 1 2 3; do
  if $SAM depth -a -Q 20 -b "$WINDOWS" "$url" \
     | awk '{k=$1"_"int($2/100000); t[k]+=$3; n[k]++}
            END{for(k in t) printf "%s\t%d\t%.4f\n", k, n[k], t[k]/n[k]}' \
     | sort > "$out.tmp" && [[ -s "$out.tmp" ]]; then
    mv "$out.tmp" "$out"; echo "  $(wc -l < "$out") windows"; exit 0
  fi
  echo "  retry $a"; rm -f "$out.tmp"; sleep $((a*30))
done
echo "FAILED $sample"; exit 1
