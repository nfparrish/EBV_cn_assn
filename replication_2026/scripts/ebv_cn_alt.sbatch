#!/usr/bin/env bash
#SBATCH --job-name=ebv_alt
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=1-300%8
#SBATCH --output=/work/nfp8/EBV/logs/alt_%A_%a.log
#
# EBV CN from an INDEPENDENT sequencing of the same lines, for test-retest against the
# NYGC 30x. Works for either the phase 3 low-coverage GRCh38DH realignment (2008-2012
# reads, different passage cells) or the 1KG-ONT hg38 alignments -- both carry chrEBV on
# the same contig, so set MANIFEST and go.
#
# DESIGN NOTE (learned the hard way): on a REMOTE CRAM, few large CONTIGUOUS regions beat
# many small scattered ones. A scattered 100 kb window still costs a seek plus a whole
# container read, so 151 scattered regions ended up pulling much of a 2-4 GB file and
# stalled against EBI. Six contiguous 2 Mb regions cost six seeks: ~60 s/sample total.
#
# Both queries are reference-free (`depth` needs POS and CIGAR, not SEQ), so no read-length
# conversion -- which matters because phase 3 mixed read lengths and sequencing centres.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
unset REF_PATH REF_CACHE

MASK_START=139700; MASK_END=151400
WT_START=144500;   WT_END=151000
MANIFEST=${MANIFEST:-$BASE/metadata/lowcov_subset.tsv}
TAG=${TAG:-lowcov}
REGIONS=(chr1:40000000-42000000 chr4:90000000-92000000 chr6:100000000-102000000
         chr11:40000000-42000000 chr14:60000000-62000000 chr18:40000000-42000000)

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$MANIFEST")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/$TAG/$sample.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/$TAG"
echo "$(date -Is) $sample ($TAG) on $(hostname)"

seek() {  # region -> "meandepth n_positions"
  local r=$1
  for a in 1 2 3; do
    if res=$($SAM depth -a -Q 20 -r "$r" "$url" 2>/dev/null \
             | awk '{t+=$3; n++} END{printf "%.5f %d", (n?t/n:0), n}') && [[ -n $res ]]; then
      echo "$res"; return 0
    fi
    sleep $((a*20))
  done
  echo "0 0"
}

# --- viral numerator, excluding the B95-8-restored interval
read -r ebv_depth ebv_n wt_depth < <($SAM depth -a -Q 20 -r chrEBV "$url" 2>/dev/null \
  | awk -v ms=$MASK_START -v me=$MASK_END -v ws=$WT_START -v we=$WT_END '
      {if ($2<ms || $2>me) {t+=$3; n++}
       if ($2>=ws && $2<=we) {w+=$3; wn++}}
      END{printf "%.5f %d %.5f\n", (n?t/n:0), n, (wn?w/wn:0)}')

# chrEBV is the LAST contig: a truncated stream silently looks like a real zero.
if [[ ${ebv_n:-0} -lt 150000 ]]; then
  echo "TRUNCATED/absent chrEBV for $sample (n=$ebv_n) -- discarding"; exit 1
fi

# --- autosomal denominator: median over the contiguous regions
vals=(); for r in "${REGIONS[@]}"; do
  read -r d n <<<"$(seek "$r")"
  [[ ${n:-0} -gt 1000000 ]] && vals+=("$d")
  echo "  $r depth=$d n=$n"
done
(( ${#vals[@]} >= 4 )) || { echo "FAILED $sample: only ${#vals[@]} usable regions"; exit 1; }
auto=$(printf '%s\n' "${vals[@]}" | sort -g | awk '{v[NR]=$1} END{print (NR%2 ? v[int(NR/2)+1] : (v[NR/2]+v[NR/2+1])/2)}')

cn=$(awk -v e="$ebv_depth" -v a="$auto" 'BEGIN{printf (a>0 ? "%.4f" : "NA"), (a>0?2*e/a:0)}')
wtf=$(awk -v w="$wt_depth" -v e="$ebv_depth" 'BEGIN{printf (e>0 ? "%.6f" : "NA"), (e>0?w/e:0)}')
printf 'sample\tdataset\tebv_depth\tauto_median_depth\tn_regions\tebv_cn\twt_fraction\n' > "$out"
printf '%s\t%s\t%.4f\t%.4f\t%d\t%s\t%s\n' "$sample" "$TAG" "$ebv_depth" "$auto" "${#vals[@]}" "$cn" "$wtf" >> "$out"
column -t "$out"
