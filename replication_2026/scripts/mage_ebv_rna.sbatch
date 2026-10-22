#!/usr/bin/env bash
#SBATCH --job-name=mage_ebv_rna
#SBATCH --account=parrishlab
#SBATCH --partition=common,scavenger
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --array=1-779%50
#SBATCH --output=/work/nfp8/EBV/logs/ebv_rna_%A_%a.log
set -euo pipefail

BASE=/work/nfp8/MAGE
OUTD=/work/nfp8/EBV/rna_mage
MM2=/hpc/group/parrishlab/envs/vcall-env/bin/minimap2
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
EBV_REF=/hpc/group/parrishlab/refs/ebv/chrEBV.fa
HOST_MMI=/work/nfp8/EBV/metadata/human_nr_tx.mmi
QUANT=/work/nfp8/EBV/scripts/ebv_rna_quant.py

mkdir -p "$OUTD/genes"
run=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$BASE/metadata/download_plan.tsv")
[[ -n "${run:-}" ]] || { echo "no row $SLURM_ARRAY_TASK_ID"; exit 1; }
[[ -f $BASE/fastq/$run/.done ]] || { echo "$run not downloaded yet"; exit 0; }
[[ -s $OUTD/${run}.ebv_rna.tsv ]] && { echo "$run already done"; exit 0; }

r1=$BASE/fastq/$run/${run}_1.fastq.gz
r2=$BASE/fastq/$run/${run}_2.fastq.gz
tmp=$(mktemp -d "${TMPDIR:-/tmp}/ebv_rna_${run}_XXXXXX")
trap 'rm -rf "$tmp"' EXIT

# 1. Align paired-end RNA-seq reads to chrEBV in splice-aware short-read mode (-ax splice:sr)
$MM2 -t "${SLURM_CPUS_PER_TASK:-4}" -ax splice:sr --secondary=no "$EBV_REF" "$r1" "$r2" 2>/dev/null \
  | $SAM view -h -F 2308 - > "$tmp/ebv.sam"

# 2. Extract candidate viral reads to FASTQ (with /1 or /2 suffix) for competitive host-transcriptome filtering
awk -F'\t' '$1 !~ /^@/ && $10 != "*" {
  suf = (and($2, 64) ? "/1" : (and($2, 128) ? "/2" : ""));
  printf "@%s%s\n%s\n+\n%s\n", $1, suf, $10, $11
}' "$tmp/ebv.sam" > "$tmp/cand.fq"

# 3. Map candidate viral reads against non-redundant human GENCODE v45 transcriptome (63k genes)
if [[ -s "$tmp/cand.fq" ]]; then
  $MM2 -t "${SLURM_CPUS_PER_TASK:-4}" -x sr --secondary=no "$HOST_MMI" "$tmp/cand.fq" > "$tmp/host.paf" 2>/dev/null
else
  : > "$tmp/host.paf"
fi

# 4. Quantify non-overlapping kinetic classes, splice junctions, and decoy filters
python3 "$QUANT" "$run" "$tmp/ebv.sam" "$tmp/host.paf" "$OUTD/${run}.ebv_rna.tsv" "$OUTD/genes/${run}.genes.tsv"
echo "Done $run"
