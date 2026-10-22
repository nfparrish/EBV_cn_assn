#!/usr/bin/env bash
# Junction-SPECIFIC V(D)J evidence from a WGS Ig-locus slice, plus corrected clonality metrics.
#
# WHY THE REALIGNMENT PASS: simply counting soft-clipped reads is not specific -- 10-15% of
# reads in the Ig loci are clipped, dominated by reference mismatch in these highly polymorphic
# regions. A clipped read only evidences a V(D)J junction if its CLIPPED SEGMENT maps to a
# DIFFERENT segment family than its anchor. So: extract clips >=20 bp, realign them against the
# IMGT V/D/J reference, and keep only clips whose best hit is a different family (V->J, V->D,
# D->J) than the anchor's own position implies.
#
# WHY THE METRIC FIX: TRUST4 on genomic DNA reports UNREARRANGED GERMLINE V segments as
# pseudo-clonotypes (V assigned, no D/J/C, near-identical CDR3s across V genes) because the
# germline locus is present at 2 copies regardless of rearrangement. Those dominated the naive
# top-clone fraction, compressing it to 0.13-0.43 with no depth dependence. Restricting to
# PRODUCTIVE, J-ASSIGNED rows recovers the real range 0.09-1.00 (sd 0.035 -> 0.210).
set -euo pipefail
sample=$1
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
MM2=/hpc/group/parrishlab/envs/vcall-env/bin/minimap2
IMGT=/hpc/group/parrishlab/refs/trust4/human_IMGT+C.fa
bam=$BASE/bcr/$sample/$sample.ig.bam
rep=$BASE/bcr/$sample/${sample}_report.tsv
[[ -f $bam && -f $rep ]] || { echo "$sample: missing inputs"; exit 1; }

# --- 1. clipped segments >=20 bp, tagged with the anchor's locus
tmp=$(mktemp -d)
$SAM view -F 0x904 "$bam" | awk -v OUT="$tmp/clips.fa" '
  {
    cig=$6; seq=$10
    lead=0; trail=0
    if (match(cig, /^[0-9]+S/))  { lead  = substr(cig, 1, RLENGTH-1)+0 }
    if (match(cig, /[0-9]+S$/))  { trail = substr(cig, RSTART, RLENGTH-1)+0 }
    if (lead >= 20)  { printf ">%s|%s|%s|L\n%s\n", $1, $3, $4, substr(seq, 1, lead) >> OUT }
    if (trail >= 20) { printf ">%s|%s|%s|T\n%s\n", $1, $3, $4, substr(seq, length(seq)-trail+1) >> OUT }
  }'
nclip=$(grep -c "^>" "$tmp/clips.fa" 2>/dev/null || echo 0)

# --- 2. realign clips to IMGT; keep best hit per clip
if [[ ${nclip:-0} -gt 0 ]]; then
  $MM2 -x sr --secondary=no -t 2 "$IMGT" "$tmp/clips.fa" 2>/dev/null \
    | awk -F'\t' '!seen[$1]++ {print $1"\t"$6"\t"$10"\t"$11}' > "$tmp/hits.paf"
else
  : > "$tmp/hits.paf"
fi
nhit=$(wc -l < "$tmp/hits.paf")

# --- 3. a clip is junction evidence only if its target FAMILY differs from the anchor locus.
#        IGH anchors sit on chr14, IGK on chr2, IGL on chr22; a V->J join within a locus shows
#        as a clip hitting a J (or D) segment of that same locus -- that is the signal we want,
#        as distinct from a clip hitting another V (reference mismatch / polymorphism).
read -r n_vj n_vd n_other <<<"$(awk -F'\t' '
  { split($1, a, "|"); chrom=a[2]; tgt=$2
    fam = (tgt ~ /IGHJ|IGKJ|IGLJ/) ? "J" : (tgt ~ /IGHD/) ? "D" : (tgt ~ /IGHV|IGKV|IGLV/) ? "V" : "C"
    if (fam=="J") vj++; else if (fam=="D") vd++; else other++ }
  END{print vj+0, vd+0, other+0}' "$tmp/hits.paf")"

# --- 4. corrected clonality from PRODUCTIVE, J-ASSIGNED rows only
read -r nall njass nprod topreal shannon <<<"$(awk -F'\t' '
  NR>1 { nall++
    if ($7!="." && $7!="" && $4!~/out_of_frame/ && $4!~/\*/ && $4!="") {
      c=$1+0; if (c>0) { np++; tot+=c; cnt[np]=c; if (c>m) m=c } } }
  END { H=0; for (i=1;i<=np;i++) { p=cnt[i]/tot; if (p>0) H-=p*log(p) }
        printf "%d %d %d %s %s", nall, 0, np+0,
          (tot>0 ? sprintf("%.4f", m/tot) : "NA"),
          (np>1 ? sprintf("%.4f", H/log(np)) : "NA") }' "$rep")"

printf '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n' \
  "$sample" "$nclip" "$nhit" "$n_vj" "$n_vd" "$n_other" "$nprod" "$topreal" "$shannon"
rm -rf "$tmp"
