#!/usr/bin/env bash
#SBATCH --job-name=ebv_gcb
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=1-2502%40
#SBATCH --output=/work/nfp8/EBV/logs/gcb_%A_%a.log
#
# Per-sample GC-coverage bias measured WITHIN chrEBV itself.
#
# Why: IR1 is 67.2% GC against 57.8% for unique chrEBV, a 9.4-point gap inside the viral genome,
# and the IR1 unit estimate correlates with the autosomal GC-bias metric (rho +0.184, p=0.0085).
# Since GC bias is itself population-structured (eta-sq 0.083), it is a candidate confound for the
# Europe-vs-rest IR1 split. An autosomal metric is only available for 199 samples and is measured
# on a different molecule; this one is measured on the SAME molecule as IR1, for all 2,502.
#
# Metric: mean depth in the top-GC quartile of unique 500 bp bins over the bottom-GC quartile
# (quartile cuts 0.616 / 0.530 from metadata/chrEBV_gc_500bp.tsv). Repeats and the B95-8 deletion
# are excluded by construction, since the bin list covers unique sequence only.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
sample=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$BASE/ebv_cn_1kGP_2504.tsv")
[[ -n "${sample:-}" ]] || exit 1
bam=$BASE/slices/$sample.chrEBV.bam
out=$BASE/gcbias/$sample.tsv
[[ -f $out ]] && exit 0
[[ -f $bam ]] || { echo "no slice for $sample"; exit 1; }
mkdir -p "$BASE/gcbias"

$SAM depth -a -Q 0 -r chrEBV "$bam" \
 | awk -v GCF="$BASE/metadata/chrEBV_gc_500bp.tsv" -v S="$sample" '
     BEGIN { while ((getline l < GCF) > 0) { split(l,a,"\t"); gc[a[1]+0]=a[3]+0 } }
     { b=int(($2-1)/500); if (!(b in gc)) next
       if (gc[b] >= 0.616) { h+=$3; nh++ } else if (gc[b] <= 0.530) { lo+=$3; nl++ }
       t+=$3; nt++ }
     END { hm=(nh?h/nh:0); lm=(nl?lo/nl:0); tm=(nt?t/nt:0)
           printf "%s\t%.4f\t%.4f\t%.4f\t%s\n", S, hm, lm, tm, (lm>0 ? sprintf("%.4f", hm/lm) : "NA") }' > "$out.tmp"
mv "$out.tmp" "$out"; cat "$out"
