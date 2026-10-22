#!/usr/bin/env bash
# Per-site mean strand bias across the cohort. Strand bias separates artifact classes
# (T>G/A>C/T>C/C>A: mean SB 13-18, 26-37% badly biased) from credible ones (G>A/A>G/indel:
# SB 3-5). It must be applied PER SITE, not per mutation class -- chrEBV:65,783 is T>C, the
# worst class, yet has SB 2.8 and co-segregates perfectly with 132,965.
set -euo pipefail
cd /work/nfp8/EBV
for f in variants/lofreq/*.vcf.gz; do zcat "$f"; done 2>/dev/null | grep -v '^#' \
 | awk -F'\t' '{ sb=0; n=split($8,I,";"); for(i=1;i<=n;i++) if (I[i] ~ /^SB=/) sb=substr(I[i],4)+0
     k=$2"\t"$4"\t"$5; s[k]+=sb; c[k]++; if (sb>=20) hi[k]++ }
   END{ for (k in c) printf "%s\t%d\t%.2f\t%.3f\n", k, c[k], s[k]/c[k], hi[k]/c[k] }' \
 | sort -k1,1n > metadata/site_strandbias.tsv
echo "wrote $(wc -l < metadata/site_strandbias.tsv) sites"
