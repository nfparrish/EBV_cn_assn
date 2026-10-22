#!/usr/bin/env bash
# Aggregate all 2,502 LoFreq VCFs into one site x sample allele-fraction table.
# With every line descended from one B95-8 ancestor, the cohort IS the panel of normals:
# a site variable at similar AF in a large fraction of samples is systematic error
# (mapping / reference), while a site variable in a few is a candidate derived variant.
set -euo pipefail
BASE=/work/nfp8/EBV
out=$BASE/variants/all_calls.tsv
: > "$out"
n=0
for f in "$BASE"/variants/lofreq/*.vcf.gz; do
  s=$(basename "$f" .vcf.gz)
  zcat "$f" | awk -v S="$s" -F'\t' '
    !/^#/ {
      af=""; dp=""
      nf=split($8, I, ";")
      for (i=1;i<=nf;i++) { if (I[i] ~ /^AF=/) af=substr(I[i],4); if (I[i] ~ /^DP=/) dp=substr(I[i],4) }
      if (af != "") print S"\t"$2"\t"$4"\t"$5"\t"af"\t"dp
    }' >> "$out"
  n=$((n+1))
done
echo "aggregated $n samples -> $(wc -l < "$out") variant records"
