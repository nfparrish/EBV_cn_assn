#!/usr/bin/env bash
#SBATCH --job-name=ebv_lowcov
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --array=1-300%8
#SBATCH --output=/work/nfp8/EBV/logs/lowcov_%A_%a.log
#
# EBV copy number from the PHASE 3 LOW-COVERAGE GRCh38DH realignment -- reads generated
# 2008-2012 and realigned 2015, versus the NYGC 30x generated 2019. Same individuals,
# DIFFERENT passage cells and different Coriell DNA shipments, so this is a genuine
# BIOLOGICAL replicate for EBV CN, which the MAGE RNA replicates cannot provide
# (CN comes from a single WGS per donor).
#
# Same GRCh38DH reference as the 30x work, so chrEBV coordinates and mapping behaviour
# are directly comparable and the test-retest is not confounded by reference choice.
#
# Uses `samtools depth` on BOTH numerator and denominator: reference-free (needs POS and
# CIGAR, not SEQ) and free of any read-length conversion, which matters here because
# phase 3 mixed read lengths and sequencing centres.
#
# Concurrency is held to 8: EBI also serves the concurrent MAGE FASTQ download.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
unset REF_PATH REF_CACHE

MASK_START=139700; MASK_END=151400
WT_START=144500;   WT_END=151000
MANIFEST=${MANIFEST:-$BASE/metadata/lowcov_subset.tsv}
WINDOWS=$BASE/metadata/windows_v2.bed

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$MANIFEST")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/lowcov/$sample.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/lowcov"
echo "$(date -Is) $sample on $(hostname)"

# one combined region file: chrEBV plus the 150 autosomal windows
reg=$(mktemp); printf 'chrEBV\t0\t171823\n' > "$reg"; cat "$WINDOWS" >> "$reg"

ok=0
for a in 1 2 3; do
  if $SAM depth -a -Q 20 -b "$reg" "$url" > "$out.depth" 2>/dev/null && [[ -s "$out.depth" ]]; then
    ok=1; break
  fi
  echo "  retry $a"; sleep $((a*30))
done
[[ $ok -eq 1 ]] || { echo "FAILED $sample"; rm -f "$reg" "$out.depth"; exit 1; }

# Guard against a truncated stream: chrEBV is the last contig, so partial output would
# silently give ebv_depth=0. Demand near-complete chrEBV and window coverage.
ebv_bp=$(awk '$1=="chrEBV"' "$out.depth" | wc -l)
win_n=$(awk '$1!="chrEBV"{print $1"_"int($2/100000)}' "$out.depth" | sort -u | wc -l)
if [[ $ebv_bp -lt 170000 || $win_n -lt 140 ]]; then
  echo "TRUNCATED $sample (chrEBV bp=$ebv_bp, windows=$win_n) -- discarding"
  rm -f "$reg" "$out.depth"; exit 1
fi

awk -v s="$sample" -v ms=$MASK_START -v me=$MASK_END -v ws=$WT_START -v we=$WT_END '
  $1=="chrEBV" {
      if ($2 < ms || $2 > me) { et += $3; en++ }
      if ($2 >= ws && $2 <= we) { wt += $3; wn++ }
      next }
  { k = $1"_"int($2/100000); at[k] += $3; an[k]++ }
  END {
      i = 0
      for (k in at) v[++i] = at[k]/an[k]
      asort(v)
      med = (i % 2 ? v[int(i/2)+1] : (v[i/2]+v[i/2+1])/2)
      ed = (en ? et/en : 0); wd = (wn ? wt/wn : 0)
      printf "%s\t%.4f\t%.4f\t%d\t%s\t%.4f\t%s\n", s, ed, med, i,
             (med>0 ? sprintf("%.4f", 2*ed/med) : "NA"), wd,
             (ed>0 ? sprintf("%.6f", wd/ed) : "NA") }' "$out.depth" > "$out"

printf 'sample\tebv_depth\tauto_median_depth\tn_windows\tebv_cn\twt_detector_depth\twt_fraction\n' > "$out.hdr"
cat "$out" >> "$out.hdr"; mv "$out.hdr" "$out"
rm -f "$reg" "$out.depth"
column -t "$out"
