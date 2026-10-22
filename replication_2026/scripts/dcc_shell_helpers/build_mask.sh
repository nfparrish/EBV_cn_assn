#!/usr/bin/env bash
# Derive an EBV mappability/repeat mask empirically from the cohort, rather than trusting
# remembered IR1-IR4 / terminal-repeat coordinates. With MAPQ>=20 filtering, multi-mapping
# repeat regions show DEPLETED coverage; collapsed repeats show inflated coverage.
# Per position, take the median across samples of (depth / that sample's own median depth).
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
cd "$BASE"
ls slices/*.chrEBV.bam | shuf -n 60 --random-source=<(yes) > /tmp/mask_bams.txt
echo "profiling $(wc -l < /tmp/mask_bams.txt) slices"
$SAM depth -a -Q 20 -r chrEBV -f /tmp/mask_bams.txt > /tmp/mask_depth.tsv
echo "positions: $(wc -l < /tmp/mask_depth.tsv)"
awk '
  NR==FNR { for (i=3; i<=NF; i++) { s[i]+=$i; c[i]++ } ; next }
  FNR==1 { for (i in s) med[i] = s[i]/c[i] }            # per-sample mean depth as the scale
  {
    n=0; delete v
    for (i=3; i<=NF; i++) if (med[i] > 5) v[++n] = $i/med[i]
    if (n < 10) next
    asort(v)
    r = (n%2 ? v[int(n/2)+1] : (v[n/2]+v[n/2+1])/2)
    printf "%d\t%.4f\n", $2, r
  }' /tmp/mask_depth.tsv /tmp/mask_depth.tsv > "$BASE/metadata/ebv_norm_depth.tsv"
echo "wrote normalised depth profile"
# contiguous runs >=200 bp outside [0.5, 2.0] become mask intervals
awk -F'\t' '
  { bad = ($2 < 0.5 || $2 > 2.0) ? 1 : 0
    if (bad && !inrun) { start=$1; inrun=1 }
    else if (!bad && inrun) { if ($1-start >= 200) print start"\t"$1; inrun=0 } }
  END { if (inrun && $1-start >= 200) print start"\t"$1 }' "$BASE/metadata/ebv_norm_depth.tsv" \
  > "$BASE/metadata/ebv_mask_raw.tsv"
echo "=== empirically derived mask intervals (>=200 bp):"
awk -F'\t' '{printf "  chrEBV:%-7d-%-7d  %6d bp\n", $1, $2, $2-$1}' "$BASE/metadata/ebv_mask_raw.tsv"
awk -F'\t' '{t+=$2-$1} END{printf "total masked: %d bp (%.1f%% of genome)\n", t, 100*t/171823}' "$BASE/metadata/ebv_mask_raw.tsv"
