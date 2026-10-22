#!/usr/bin/env bash
#SBATCH --job-name=bcr_ont
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=1-907%6
#SBATCH --output=/work/nfp8/EBV/logs/bcront_%A_%a.log
#
# BCR rearrangements from ONT long reads, for comparison against the 30x short-read calls.
#
# WHY NOT TRUST4 HERE: TRUST4 assembles short reads. With long reads no assembly is needed --
# a single read spans the whole V(D)J, so each read IS a rearrangement call. Assembly would
# only add error.
#
# WHY SEGMENT USAGE RATHER THAN EXACT CDR3: this ONT data is R9.4.1 basecalled with Guppy 6
# (~95-97% accurate), so a ~60 bp CDR3 carries 2-3 errors. Exact CDR3 nucleotide matching
# across platforms would fail for technical reasons alone. V/J SEGMENT IDENTITY is robust to
# that error rate, so the cross-platform comparison is built on V-J pair usage.
#
# A read is rearrangement evidence when it aligns to BOTH a V and a J segment of the same
# locus -- which germline-configuration DNA cannot produce, since V and J are megabases apart
# until recombination joins them.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
MM2=/hpc/group/parrishlab/envs/vcall-env/bin/minimap2
IMGT=/hpc/group/parrishlab/refs/trust4/human_IMGT+C.fa
# CRAM->FASTA needs SEQ reconstruction, so the human reference must be supplied with -T
HREF=/hpc/group/parrishlab/refs/GRCh38_full/GRCh38_full_analysis_set_plus_decoy_hla.fa
PAD=50000
unset REF_PATH REF_CACHE

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$BASE/metadata/ont_overlap.tsv")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/bcr_ont/$sample.vj.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/bcr_ont"
echo "$(date -Is) $sample on $(hostname)"

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for a in 1 2 3; do
  if $SAM view -b -T "$HREF" -o "$tmp/ig.bam" "$url" \
        chr14:$((105583731-PAD))-$((106879812+PAD)) \
        chr2:$((88857161-PAD))-$((90235370+PAD)) \
        chr22:$((22026076-PAD))-$((22923034+PAD)); then
    $SAM fasta "$tmp/ig.bam" > "$tmp/reads.fa" 2>/dev/null
    [[ -s $tmp/reads.fa ]] && break
  fi
  echo "  retry $a"; sleep $((a*30))
done
[[ -s $tmp/reads.fa ]] || { echo "FAILED to slice $sample"; exit 1; }
nreads=$(grep -c "^>" "$tmp/reads.fa")
echo "  $nreads reads over the Ig loci"

# All alignments, not just the best: one read may hit V, D and J.
# NB: do NOT filter on MAPQ here. Asking for secondaries (-N 20 -p 0.5) means multi-mapping
# alignments get MAPQ 0 -- 91% of hits -- so a MAPQ filter silently discards exactly the
# multi-segment evidence being sought. Filter on matching bases and block length instead.
$MM2 -x map-ont -t "${SLURM_CPUS_PER_TASK:-4}" -N 20 -p 0.5 "$IMGT" "$tmp/reads.fa" 2>/dev/null \
  | awk -F'\t' '$10>=50 && $11>=50 {print $1"\t"$6}' | sort -u > "$tmp/hits.tsv"

# a read hitting both a V and a J of the same locus evidences a rearrangement
awk -F'\t' -v S="$sample" '
  { read=$1; seg=$2
    locus = (seg ~ /^IGH/) ? "IGH" : (seg ~ /^IGK/) ? "IGK" : (seg ~ /^IGL/) ? "IGL" : "other"
    if (locus=="other") next
    if (seg ~ /V[0-9]/) { v[read"|"locus] = seg }
    else if (seg ~ /J[0-9]/) { j[read"|"locus] = seg } }
  END {
    for (k in v) if (k in j) { split(k, a, "|"); pair[a[2]"\t"v[k]"\t"j[k]]++ }
    for (p in pair) printf "%s\t%s\t%d\n", S, p, pair[p] }' "$tmp/hits.tsv" | sort -k5,5nr > "$out.tmp"

mv "$out.tmp" "$out"
n=$(wc -l < "$out")
tot=$(awk -F'\t' '{t+=$5} END{print t+0}' "$out")
top=$(awk -F'\t' -v T="$tot" 'NR==1 && T>0 {printf "%.3f", $5/T}' "$out")
echo "  V-J pairs=$n  supporting reads=$tot  top-pair fraction=${top:-NA}"
printf '%s\t%d\t%d\t%d\t%s\n' "$sample" "$nreads" "$n" "$tot" "${top:-NA}" >> "$BASE/bcr_ont_summary.tsv"
