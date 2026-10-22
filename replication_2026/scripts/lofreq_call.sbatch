#!/usr/bin/env bash
#SBATCH --job-name=ebv_lofreq
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-2502%25
#SBATCH --output=/work/nfp8/EBV/logs/lofreq_%A_%a.log
#
# Intrahost EBV variant calling with LoFreq (best-benchmarked tool for viral intrahost
# populations; Mutect2 mitochondria-mode runs later as the orthogonal check).
#
# TWO DESIGN DECISIONS THAT MATTER MORE THAN THE CALLER:
#
# 1. DOWNSAMPLE TO A COMMON DEPTH. Viral depth spans 51x-10,892x (200-fold) and CN
#    correlates with population (eta-sq 0.126). Published benchmarks report LoFreq and
#    Mutect2 are strongly affected by copy number, so CN-dependent detection power would
#    manufacture population-structured variants -- mimicking exactly the provenance signal
#    we are testing for. Every sample is downsampled to ~150x, which retains 2,479/2,502.
#
# 2. MASK THE REPEATS, derived empirically from cohort coverage rather than from
#    remembered coordinates. The mask recovers IR1 (12.3-35.0 kb), IR2, IR3 (the EBNA1
#    Gly-Ala repeat), the terminal repeats, and re-derives the B95-8 deletion to within
#    24 bp of the independent estimate. 22.3% of the genome is excluded.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
LOFREQ=/hpc/group/parrishlab/envs/vcall-env/bin/lofreq
BGZIP=/hpc/group/parrishlab/envs/samtools-env/bin/bgzip
BCFTOOLS=/hpc/group/parrishlab/envs/samtools-env/bin/bcftools
REF=/hpc/group/parrishlab/refs/ebv/chrEBV.fa
TARGET=150

sample=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$BASE/ebv_cn_1kGP_2504.tsv")
depth=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $2}' "$BASE/ebv_cn_1kGP_2504.tsv")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
bam=$BASE/slices/$sample.chrEBV.bam
out=$BASE/variants/lofreq/$sample.vcf.gz
[[ -f $out ]] && { echo "$sample done"; exit 0; }
[[ -f $bam ]] || { echo "no slice for $sample"; exit 1; }
mkdir -p "$BASE/variants/lofreq" "$BASE/variants/tmp"
echo "$(date -Is) $sample (chrEBV depth ${depth}x) on $(hostname)"

work=$BASE/variants/tmp/$sample
# Downsample to the common target so detection power does not track CN.
frac=$(awk -v d="$depth" -v t=$TARGET 'BEGIN{printf "%.5f", (d>t ? t/d : 1)}')
if awk -v f="$frac" 'BEGIN{exit !(f<1)}'; then
  echo "  downsampling to ~${TARGET}x (fraction $frac)"
  $SAM view -b -s "$frac" -o "$work.ds.bam" "$bam"
else
  echo "  depth below target; using all reads (flagged)"
  cp "$bam" "$work.ds.bam"
fi
$SAM index "$work.ds.bam"
actual=$($SAM depth -a -Q 20 -r chrEBV "$work.ds.bam" | awk '{t+=$3;n++} END{printf "%.1f",(n?t/n:0)}')

# indel qualities, then call SNVs+indels over callable regions only
$LOFREQ indelqual --dindel -f "$REF" -o "$work.iq.bam" "$work.ds.bam" 2>/dev/null
$SAM index "$work.iq.bam"
$LOFREQ call --call-indels -f "$REF" -l "$BASE/metadata/ebv_callable.bed" \
    -q 20 -Q 20 -m 20 -o "$work.vcf" "$work.iq.bam" 2>/dev/null

$BGZIP -f -c "$work.vcf" > "$out"; $BCFTOOLS index -f "$out"
n=$(zcat "$out" | grep -vc '^#' || true)
echo "  effective depth ${actual}x, variants called: $n"
printf '%s\t%s\t%s\t%s\n' "$sample" "$depth" "$actual" "$n" >> "$BASE/variants/lofreq_summary.tsv"
rm -f "$work".*
