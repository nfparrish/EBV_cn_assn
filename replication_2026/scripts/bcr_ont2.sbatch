#!/usr/bin/env bash
#SBATCH --job-name=bcr_ont
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=1-907%6
#SBATCH --output=/work/nfp8/EBV/logs/bcront_%A_%a.log
#
# V(D)J rearrangements from ONT long reads, read off the GENOMIC alignment structure.
#
# Why this rather than aligning reads to IMGT segments: IGHJ segments are only ~50 bp and this
# data is R9.4.1 (~5% error), so short-segment alignment is unreliable -- it returned 1-read
# "clones", two of them spurious hits to chr16 IGHV orphons. Alignment STRUCTURE is robust to
# base-level error.
#
# The signature: recombination joins a V to a J ~70-900 kb away in the reference and deletes the
# interval, so a read crossing the junction aligns in two pieces -- primary + SUPPLEMENTARY --
# one in the IGHV array, one at IGHD/IGHJ. GRCh38 chr14 coordinates:
#   IGHD 105,836,765-105,845,678   IGHJ 105,863,198-105,865,678   IGHV 105,939,756-106,879,812
#
# TWO FILTERS THAT MATTER:
#   - EXCLUDE secondary alignments (0x100). The IGHV array is full of near-identical genes, so
#     secondaries generate spurious V-to-V "pairs". Keep primary + supplementary only.
#   - REQUIRE DISJOINT QUERY INTERVALS: the two alignments must cover different parts of the
#     read, which is what a genuine junction-spanning read looks like.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
HREF=/hpc/group/parrishlab/refs/GRCh38_full/GRCh38_full_analysis_set_plus_decoy_hla.fa
unset REF_PATH REF_CACHE
JD_LO=105836765; JD_HI=105865678; V_LO=105939756; V_HI=106879812

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$BASE/metadata/ont_overlap.tsv")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/bcr_ont/$sample.vdj.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/bcr_ont"
echo "$(date -Is) $sample on $(hostname)"

tmp=$(mktemp -d -p "$BASE/scratch"); trap 'rm -rf "$tmp"' EXIT
for a in 1 2 3; do
  $SAM view -b -T "$HREF" -o "$tmp/igh.bam" "$url" chr14:${JD_LO}-${V_HI} && break
  echo "  retry $a"; sleep $((a*30))
done
[[ -s $tmp/igh.bam ]] || { echo "FAILED to slice $sample"; exit 1; }

# primary + supplementary only; record query-interval start/end per alignment
$SAM view -F 0x100 "$tmp/igh.bam" | awk -v OFS='\t' '
  { cig=$6; qs=0; qlen=0
    if (match(cig,/^[0-9]+[SH]/)) qs = substr(cig,1,RLENGTH-1)+0
    c=cig; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
      l=substr(c,1,RLENGTH-1)+0; o=substr(c,RLENGTH,1); c=substr(c,RLENGTH+1)
      if (o=="M"||o=="I"||o=="="||o=="X") qlen+=l }
    print $1, $4, qs, qs+qlen }' | sort -k1,1 -k2,2n > "$tmp/aln.tsv"

awk -v OFS='\t' -v S="$sample" -v jl=$JD_LO -v jh=$JD_HI -v vl=$V_LO -v vh=$V_HI '
  { r[$1] = r[$1] $2","$3","$4";" }
  END {
    for (read in r) {
      n = split(r[read], A, ";")
      for (i=1; i<n; i++) for (j=i+1; j<n; j++) {
        split(A[i], a, ","); split(A[j], b, ",")
        p1=a[1]+0; s1=a[2]+0; e1=a[3]+0; p2=b[1]+0; s2=b[2]+0; e2=b[3]+0
        if (e1 <= s2 + 50 || e2 <= s1 + 50) {            # disjoint query intervals
          lo = (p1<p2?p1:p2); hi = (p1<p2?p2:p1)
          if (lo>=jl && lo<=jh && hi>=vl && hi<=vh) hits[lo"\t"hi]++
        } } }
    for (k in hits) print S, k, hits[k] }' "$tmp/aln.tsv" | sort -k4,4nr > "$out.tmp"
mv "$out.tmp" "$out"

n=$(wc -l < "$out"); tot=$(awk -F'\t' '{t+=$4} END{print t+0}' "$out")
top=$(awk -F'\t' -v T="$tot" 'NR==1 && T>0 {printf "%.3f", $4/T}' "$out")
echo "  junction-spanning reads=$tot across $n distinct V->J breakpoints; top=${top:-NA}"
printf '%s\t%d\t%d\t%s\n' "$sample" "$n" "$tot" "${top:-NA}" >> "$BASE/bcr_ont_summary.tsv"
