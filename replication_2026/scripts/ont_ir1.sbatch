#!/usr/bin/env bash
#SBATCH --job-name=ont_ir1
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --array=1-907%6
#SBATCH --output=/work/nfp8/EBV/logs/ontir1_%A_%a.log
#
# IR1 unit count and sequence from ONT long reads -- the orthogonal test of the short-read
# depth-ratio estimate, and the only route to a positive control for variant calling.
#
# WHY IT IS DECISIVE: short-read depth gives IR1 = 9.9-10.5 units in all four European
# populations vs 12.0-13.2 everywhere else (eta-sq 0.36, fully explained by continental group).
# Two explanations survive: a two-stock provenance model, or a technical artifact in
# repeat-region coverage correlated with sample batch. ONT has different chemistry and NO
# amplification, and counts units directly from read structure rather than from depth, so it
# cannot share the artifact. If ONT reproduces the Europe-vs-rest split, the split is real.
#
# WHY A POSITIVE CONTROL IS NEEDED: the only published B95-8 minor variants (J Virol 2017,
# PMC5686732, NC_007605 coordinates) sit at 12,540 / 12,647 / 13,210 / 13,290 / 14,450 -- all
# inside IR1, where our MAPQ>=20 SNV pipeline has 0.000-0.135 normalised depth. We currently
# cannot check a single call against external truth. ONT reads span 3-6 W units (3,072 bp each),
# which is how the 2017 authors resolved IR1 with PacBio.
#
# Two independent unit estimates per sample:
#   (a) depth ratio at MAPQ 0 over IR1 vs unique chrEBV -- same statistic as the short reads;
#   (b) per-read aligned span across IR1 / 3,072 -- structural, depth-free, ONT-only.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
HREF=/hpc/group/parrishlab/refs/GRCh38_full/GRCh38_full_analysis_set_plus_decoy_hla.fa
IR1_LO=12301; IR1_HI=35037; UNIT=3072
unset REF_PATH REF_CACHE

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "$BASE/metadata/ont_overlap.tsv")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
out=$BASE/ont_ir1/$sample.tsv
[[ -f $out ]] && { echo "$sample done"; exit 0; }
mkdir -p "$BASE/ont_ir1" "$BASE/ont_ir1/bam"
echo "$(date -Is) $sample on $(hostname)"

bam=$BASE/ont_ir1/bam/$sample.chrEBV.bam
if [[ ! -f $bam ]]; then
  for a in 1 2 3; do
    $SAM view -b -T "$HREF" -o "$bam.tmp" "$url" chrEBV && { mv "$bam.tmp" "$bam"; break; }
    echo "  retry $a"; rm -f "$bam.tmp"; sleep $((a*30))
  done
fi
[[ -f $bam ]] || { echo "FAILED to slice $sample"; exit 1; }
$SAM index "$bam"

# (a) depth-ratio estimate, MAPQ 0, directly comparable to the short-read number
read -r ir1d unid <<<"$($SAM depth -a -Q 0 -r chrEBV "$bam" | awk -v il=$IR1_LO -v ih=$IR1_HI '
  { p=$2; d=$3
    if (p>=il && p<=ih) { i+=d; ni++ }
    if ((p<12301)||(p>35037&&p<38382)||(p>39733&&p<95984)||(p>96575&&p<139724)||(p>151551&&p<169740)) { u+=d; nu++ } }
  END{ printf "%.4f %.4f\n", (ni?i/ni:0), (nu?u/nu:0) }')"
ratio=$(awk -v a="$ir1d" -v b="$unid" 'BEGIN{printf (b>0?"%.4f":"NA"), (b>0?a/b*(22737/3072):0)}')

# (b) structural estimate: reads whose alignment spans IR1, span / unit length
read -r nspan medunits maxunits <<<"$($SAM view -F 0x904 "$bam" | awk -v il=$IR1_LO -v ih=$IR1_HI -v U=$UNIT '
  { pos=$4; cig=$6; reflen=0
    while (match(cig, /^[0-9]+[MIDNSHP=X]/)) { l=substr(cig,1,RLENGTH-1)+0; o=substr(cig,RLENGTH,1); cig=substr(cig,RLENGTH+1)
      if (o=="M"||o=="D"||o=="N"||o=="="||o=="X") reflen+=l }
    s=pos; e=pos+reflen
    ov = (e<ih?e:ih) - (s>il?s:il)
    if (ov >= 2*U) { n++; v[n]=ov/U } }
  END{ if (!n) { print 0, "NA", "NA"; exit }
       asort(v); m=(n%2?v[int(n/2)+1]:(v[n/2]+v[n/2+1])/2)
       printf "%d %.2f %.2f\n", n, m, v[n] }')"

printf 'sample\tir1_depth\tunique_depth\tunits_by_depth\tn_spanning_reads\tmedian_units_spanned\tmax_units_spanned\n' > "$out"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$ir1d" "$unid" "$ratio" "$nspan" "$medunits" "$maxunits" >> "$out"
tail -1 "$out"
rm -f "$bam" "$bam.bai"     # slices are large; metrics are what we keep
