#!/usr/bin/env bash
#SBATCH --job-name=bcr_wgs
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=1-731%12
#SBATCH --output=/work/nfp8/EBV/logs/bcr_%A_%a.log
#
# BCR clonality from short-read WGS: localise the three immunoglobulin loci from the NYGC
# 30x CRAMs and extract V(D)J evidence.
#
# WHY: LCLs start polyclonal and become oligoclonal, so the repertoire is a per-line
# culture-history probe -- the covariate IGSR and Coriell do not publish. The DNA clones
# (2019 WGS) compared against the RNA clones (MAGE culture) give a per-sample measure of
# how far the cultures behind our predictor (CN) and our outcome (lytic score) have drifted,
# which bounds H1 independent of biology.
#
# WHY NOT DEPTH: an earlier idea was to read clonality off IGH coverage (the V-to-J interval
# is deleted in a clone). That fails -- LCLs are rarely monoclonal, co-resident clones delete
# DIFFERENT intervals, and their depth effects superimpose into a smear; IGH is also among the
# worst-mapping, most structurally polymorphic loci in the genome. Junction evidence instead.
#
# Loci derived from TRUST4's own segment coordinates (not hardcoded), padded 50 kb:
#   IGH chr14:105,583,731-106,879,812   IGK chr2:88,857,161-90,235,370
#   IGL chr22:22,026,076-22,923,034
#
# TRUST4 runs on the slice via `-b`, the SAME tool and IMGT reference used on the MAGE RNA,
# so DNA and RNA clonotypes share CDR3 definitions and are directly comparable.
#
# SENSITIVITY, stated honestly: a clone at frequency f yields ~30*f junction-spanning reads,
# so f >~ 0.10 is confident and ~0.05 marginal. This finds DOMINANT clones -- sufficient for
# "is it the same clone across assays" -- and will not reconstruct a polyclonal repertoire.
set -euo pipefail
BASE=/work/nfp8/EBV
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
T4=/hpc/group/parrishlab/envs/vcall-env/bin/run-trust4
HREF=/hpc/group/parrishlab/refs/GRCh38_full/GRCh38_full_analysis_set_plus_decoy_hla.fa
T4REF=/hpc/group/parrishlab/refs/trust4
PAD=50000

read -r sample url < <(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n{print $1, $2}' "${MANIFEST:-$BASE/metadata/bcr_manifest.tsv}")
[[ -n "${sample:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
outd=$BASE/bcr/$sample
[[ -f $outd/${sample}_report.tsv ]] && { echo "$sample done"; exit 0; }
mkdir -p "$outd"
echo "$(date -Is) $sample on $(hostname)"

bam=$outd/$sample.ig.bam
if [[ ! -f $bam ]]; then
  for a in 1 2 3; do
    if $SAM view -b -T "$HREF" -o "$bam.tmp" "$url" \
         chr14:$((105583731-PAD))-$((106879812+PAD)) \
         chr2:$((88857161-PAD))-$((90235370+PAD)) \
         chr22:$((22026076-PAD))-$((22923034+PAD)); then
      mv "$bam.tmp" "$bam"; break
    fi
    echo "  retry $a"; rm -f "$bam.tmp"; sleep $((a*30))
  done
fi
[[ -f $bam ]] || { echo "FAILED to slice $sample"; exit 1; }
$SAM index "$bam"
echo "  sliced $(du -h "$bam" | cut -f1), $($SAM view -c "$bam") reads"

# --- raw V(D)J evidence counts. More robust than assembly at 30x, so recorded separately.
#     softclip: >=20 bp clipped (candidate junction-spanning read)
#     discord : mate-pair spanning an apparent deletion (V-to-J excision), >1 kb insert
#     split   : supplementary alignment (SA tag)
read -r nsc ndisc nsplit <<<"$($SAM view -F 0x904 "$bam" | awk '
  {
    if ($6 ~ /^[0-9]+S/ || $6 ~ /[0-9]+S$/) {
      n=0; if (match($6,/^([0-9]+)S/,m)) n=m[1]; if (match($6,/([0-9]+)S$/,m2) && m2[1]>n) n=m2[1]
      if (n >= 20) sc++
    }
    if ($7=="=" && ($9>1000 || $9<-1000)) disc++
    for (i=12;i<=NF;i++) if ($i ~ /^SA:Z:/) { sp++; break }
  }
  END{print sc+0, disc+0, sp+0}')"
echo "  evidence: softclip>=20bp=$nsc  discordant(>1kb)=$ndisc  split(SA)=$nsplit"

# --- clonotypes, same tool/reference as the RNA side
cd "$outd"
$T4 -b "$bam" -f "$T4REF/hg38_bcrtcr.fa" --ref "$T4REF/human_IMGT+C.fa" \
    -o "$sample" -t "${SLURM_CPUS_PER_TASK:-4}" --abnormalUnmapFlag 2>&1 | tail -3

rep=$outd/${sample}_report.tsv
if [[ -f $rep ]]; then
  nclone=$(tail -n +2 "$rep" | wc -l)
  igh=$(awk -F'\t' 'NR>1 && $0 ~ /IGH/ {c++} END{print c+0}' "$rep")
  top=$(awk -F'\t' 'NR>1 {s+=$1; if ($1>m) m=$1} END{printf "%.3f", (s>0 ? m/s : 0)}' "$rep")
  echo "  clonotypes=$nclone IGH=$igh top-clone=$top"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$nsc" "$ndisc" "$nsplit" "$nclone" "$igh" "$top" \
    >> "$BASE/bcr_wgs_summary.tsv"
else
  echo "  WARNING: no TRUST4 report for $sample"
fi
