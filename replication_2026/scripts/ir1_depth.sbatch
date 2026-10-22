#!/usr/bin/env bash
#SBATCH --job-name=ebv_ir1
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=1-2502%40
#SBATCH --output=/work/nfp8/EBV/logs/ir1_%A_%a.log
#
# IR1 (BamHI-W major internal repeat) copy number per line, WITHOUT the MAPQ filter.
#
# Rationale: IR1 is a tandem array collapsed to a single copy in NC_007605, so reads from all
# N units pile onto that one copy. At MAPQ>=20 they are discarded as multi-mapping and IR1 reads
# ~0x depth -- correct for SNV calling, wrong here, because for repeat-unit counting the
# multi-mapping reads ARE the signal. Depth(IR1)/Depth(unique) therefore estimates unit count.
#
# Why it matters: J Virol 2017 (PMC5686732, NC_007605 coordinates -- same reference as ours)
# reports IR1 unit count differing BETWEEN B95-8 preparations: 11.6 in the B95-8 cell line vs
# 6.6 in the HB9 BAC clone. IR1 copy number is therefore a documented stock-to-stock difference,
# which makes it an independent provenance marker in exactly the region where stock differences
# are known -- and the one region our SNV pipeline masks out.
#
# No network: operates on the chrEBV slices already on disk, which were extracted WITHOUT a MAPQ
# filter, so the multi-mapping reads are present.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
IR1_LO=12301; IR1_HI=35037          # empirical IR1 bounds
CP_LO=10000;  CP_HI=11500           # Cp promoter region, unique, immediately 5' of IR1

sample=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$BASE/ebv_cn_1kGP_2504.tsv")
[[ -n "${sample:-}" ]] || exit 1
bam=$BASE/slices/$sample.chrEBV.bam
out=$BASE/ir1/$sample.tsv
[[ -f $out ]] && exit 0
[[ -f $bam ]] || { echo "no slice for $sample"; exit 1; }
mkdir -p "$BASE/ir1"

# MAPQ 0 (-Q 0): keep multi-mapping reads. Unique denominator avoids all repeat regions.
read -r ir1 cp uniq <<<"$($SAM depth -a -Q 0 -r chrEBV "$bam" | awk -v il=$IR1_LO -v ih=$IR1_HI \
  -v cl=$CP_LO -v ch=$CP_HI '
   { p=$2; d=$3
     if (p>=il && p<=ih) { i+=d; ni++ }
     else if (p>=cl && p<=ch) { c+=d; nc++ }
     # unique sequence: outside every repeat and outside the B95-8 deletion
     if ((p<12301) || (p>35037 && p<38382) || (p>39733 && p<95984) || (p>96575 && p<139724) \
         || (p>151551 && p<169740)) { u+=d; nu++ } }
   END{ printf "%.4f %.4f %.4f\n", (ni?i/ni:0), (nc?c/nc:0), (nu?u/nu:0) }')"

ratio=$(awk -v a="$ir1" -v b="$uniq" 'BEGIN{printf (b>0 ? "%.3f" : "NA"), (b>0?a/b:0)}')
cpr=$(awk -v a="$ir1" -v b="$cp" 'BEGIN{printf (b>0 ? "%.3f" : "NA"), (b>0?a/b:0)}')
printf 'sample\tir1_depth_mapq0\tcp_depth\tunique_depth\tir1_copies_vs_unique\tir1_copies_vs_Cp\n' > "$out"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$ir1" "$cp" "$uniq" "$ratio" "$cpr" >> "$out"
cat "$out" | tail -1
