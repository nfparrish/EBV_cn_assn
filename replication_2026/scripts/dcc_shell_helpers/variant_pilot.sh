#!/usr/bin/env bash
# Pilot the provenance test: joint-call EBV variants across IBS and GBR lines.
# These two populations share European ancestry but differ ~4x in EBV copy number
# (IBS ~100 vs GBR ~24 copies/cell), so they are the sharpest available contrast
# between "LCL provenance" and "ancestry".
set -euo pipefail
BASE=/work/nfp8/EBV
BIN=/hpc/group/parrishlab/envs/samtools-env/bin
REF=/hpc/group/parrishlab/refs/ebv/chrEBV.fa
META=/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt
mkdir -p "$BASE/variants"

for pop in IBS GBR; do
  awk -F'\t' -v p="$pop" 'NR>1 && $7==p {print $4}' "$META" | sort -u \
    | while read -r s; do [[ -f $BASE/slices/$s.chrEBV.bam ]] && echo "$s"; done > "$BASE/variants/$pop.ids"
  echo "$pop: $(wc -l < "$BASE/variants/$pop.ids") lines with slices"
done

cat "$BASE/variants/IBS.ids" "$BASE/variants/GBR.ids" > "$BASE/variants/pilot.ids"
sed 's|^|'"$BASE"'/slices/|; s|$|.chrEBV.bam|' "$BASE/variants/pilot.ids" > "$BASE/variants/pilot.bams"

# Exclude the B95-8-restored interval (no template) and the ~143.5-144 kb mismapping blip.
printf 'chrEBV\t1\t139699\nchrEBV\t151401\t171823\n' > "$BASE/variants/callable.bed"

echo "joint-calling $(wc -l < "$BASE/variants/pilot.bams") samples..."
$BIN/bcftools mpileup -f "$REF" -b "$BASE/variants/pilot.bams" -R "$BASE/variants/callable.bed" \
    -a AD,DP -q 20 -Q 20 -d 2000 -Ou 2>/dev/null \
  | $BIN/bcftools call -mv -Ou 2>/dev/null \
  | $BIN/bcftools norm -f "$REF" -Ou 2>/dev/null \
  | $BIN/bcftools filter -i 'QUAL>=50 && INFO/DP>=1000' -Oz -o "$BASE/variants/pilot.vcf.gz" 2>/dev/null
$BIN/bcftools index -f "$BASE/variants/pilot.vcf.gz"
echo "variant sites called: $($BIN/bcftools view -H "$BASE/variants/pilot.vcf.gz" | wc -l)"
