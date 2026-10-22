#!/usr/bin/env bash
#SBATCH --job-name=ebv_slice
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --array=1-2504%15
#SBATCH --output=/work/nfp8/EBV/logs/%A_%a.log
#
# Per 1kGP 30x sample: remote-slice chrEBV from the NYGC CRAM over HTTPS (never
# downloading the whole 15.4 GB file), keep the viral reads for later variant calling,
# and emit two EBV copy-number estimates plus a B95-8-purity metric.
#
# chrEBV == NC_007605/AJ507799.2: B95-8 with the ~12 kb deletion restored from Raji.
# Measured empirically (NA12718, 4 kb then 500 bp bins) that restored interval is
# chrEBV:139,700-151,400 and reads ~0x in a pure-B95-8 line, hence:
#   MASK  139700-151400  excluded from the CN numerator (no template in B95-8)
#   WTDET 144500-151000  wild-type co-infection detector; clean-zero sub-interval,
#                        chosen to skip a mismapping blip at ~143.5-144 kb
#
# No human reference is needed: the viral slice decodes against a 172 kb chrEBV-only
# FASTA, and the autosomal denominator is a reference-free record count. (Relying on
# htslib's REF_PATH registry at EBI was the original bottleneck -- ~10 min/sample.)
#
# Resumable and idempotent: samples with a cn/.done marker are skipped, so the array
# can be resubmitted to retry failures. Override MANIFEST for a different cohort.
set -euo pipefail

BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
EBVREF=/hpc/group/parrishlab/refs/ebv/chrEBV.fa
unset REF_PATH REF_CACHE

EBV_LEN=171823
MASK_START=139700; MASK_END=151400
WT_START=144500;   WT_END=151000
MAPQ=20
MANIFEST=${MANIFEST:-$BASE/metadata/cram_manifest.tsv}
WINDOWS=$BASE/metadata/autosomal_windows.bed
WINDOW_BP=$(awk '{t+=$3-$2} END{print t}' "$WINDOWS")

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$MANIFEST")
[[ -n "${sample:-}" ]] || { echo "no manifest row for task $SLURM_ARRAY_TASK_ID"; exit 1; }
url=${url/ftp:\/\//https://}

[[ -f $BASE/cn/.done/$sample ]] && { echo "$sample already done"; exit 0; }
mkdir -p "$BASE/slices" "$BASE/cn/.done"
echo "$(date -Is) $sample on $(hostname)"

retry() { for a in 1 2 3 4 5; do "$@" && return 0; echo "  retry $a"; sleep $((a*30)); done; return 1; }

vbam=$BASE/slices/$sample.chrEBV.bam

echo "  [$(date +%T)] viral slice"
retry $SAM view -b -T "$EBVREF" -o "$vbam.tmp" "$url" chrEBV
mv "$vbam.tmp" "$vbam"; $SAM index "$vbam"

echo "  [$(date +%T)] autosomal read count"
acount=$BASE/cn/$sample.autocount
retry bash -c "$SAM view -c -q $MAPQ -M -L '$WINDOWS' '$url' > '$acount'"
auto_reads=$(cat "$acount"); rm -f "$acount"
[[ -n "$auto_reads" && "$auto_reads" -gt 0 ]] || { echo "autosomal count failed"; exit 1; }

echo "  [$(date +%T)] viral depths"
# Mean chrEBV depth excluding the restored interval, and the wild-type detector depth.
ebv_depth=$($SAM depth -a -Q $MAPQ -r chrEBV "$vbam" \
  | awk -v s=$MASK_START -v e=$MASK_END '$2<s || $2>e {t+=$3; n++} END{print (n?t/n:0)}')
wt_depth=$($SAM depth -a -Q $MAPQ -r "chrEBV:$WT_START-$WT_END" "$vbam" \
  | awk '{t+=$3; n++} END{print (n?t/n:0)}')
# Reads outside the mask, for the independent read-count estimator
ebv_reads=$($SAM view -c -q $MAPQ "$vbam" "chrEBV:1-$((MASK_START-1))" "chrEBV:$((MASK_END+1))-$EBV_LEN")
read_len=$($SAM view "$vbam" | awk '{t+=length($10); n++} END{print (n?t/n:0)}')

eff_bp=$((EBV_LEN - (MASK_END - MASK_START)))
read -r cn_depth cn_reads auto_depth wt_frac < <(awk -v ed="$ebv_depth" -v er="$ebv_reads" \
  -v ar="$auto_reads" -v rl="$read_len" -v eb="$eff_bp" -v wb="$WINDOW_BP" -v wd="$wt_depth" 'BEGIN{
    ad = ar * rl / wb;                       # autosomal mean depth
    printf "%.3f %.3f %.2f %.5f\n",
      (ad>0 ? 2*ed/ad : -1),                 # depth-ratio CN
      (ar>0 ? 2*(er/eb)/(ar/wb) : -1),       # read-count-ratio CN (read length cancels)
      ad,
      (ed>0 ? wd/ed : -1)}')

res=$BASE/cn/$sample.tsv
printf 'sample\tebv_depth\tauto_depth\tebv_cn_depth_ratio\tebv_cn_read_ratio\tebv_reads\tauto_reads\tread_len\twt_detector_depth\twt_fraction\n' > "$res"
printf '%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.3f\t%s\n' \
  "$sample" "$ebv_depth" "$auto_depth" "$cn_depth" "$cn_reads" \
  "$ebv_reads" "$auto_reads" "$read_len" "$wt_depth" "$wt_frac" >> "$res"
column -t "$res"
touch "$BASE/cn/.done/$sample"
echo "  [$(date +%T)] done"
