#!/usr/bin/env bash
# Is the RNA-junction-into-ONT search a method failure or a POWER failure?
# Positive control: a germline V segment must be present in every B cell, so it must be found.
set -euo pipefail
cd /work/nfp8/EBV/bcr_blast/proto
BL=/hpc/group/parrishlab/envs/vcall-env/bin
IMGT=/hpc/group/parrishlab/refs/trust4/human_IMGT+C.fa

# positive control: pull one germline IGHV and one IGKV from the IMGT reference
python3 -I - <<'PY'
want=["IGHV3-23*01","IGKV3-20*01"]; keep=False; out=open("poscontrol.fa","w"); cur=None
for line in open("/hpc/group/parrishlab/refs/trust4/human_IMGT+C.fa"):
    if line.startswith(">"):
        cur=line[1:].split()[0]; keep = cur in want
        if keep: out.write(f">{cur}_germline\n")
    elif keep: out.write(line)
out.close()
PY
echo "positive-control queries: $(grep -c '^>' poscontrol.fa)"

echo
printf '%-34s %10s %10s\n' "query set / database" "alignments" "queries_hit"
for q in poscontrol.fa HG00096_junction_queries.fa; do
  nq=$(grep -c '^>' "$q")
  for d in HG00096 HG00130; do
    tag=$([ "$d" = HG00096 ] && echo OWN || echo CTRL)
    tmp=$(mktemp)
    $BL/blastn -query "$q" -db ${d}_db -outfmt '6 qseqid pident length qlen' \
       -evalue 1e-5 -word_size 11 -num_threads 2 2>/dev/null \
       | awk '$2>=85 && $3>=0.6*$4' > "$tmp" || true
    printf '%-34s %10s %10s\n' "$(basename "$q" .fa) -> $d($tag)" "$(wc -l < "$tmp")" "$(cut -f1 "$tmp" | sort -u | wc -l)/$nq"
    rm -f "$tmp"
  done
done

echo
echo "=== ONT IG coverage, which sets the power ceiling"
for d in HG00096 HG00130; do
  python3 -I - "$d" <<'PY'
import sys
d=sys.argv[1]; n=0; tot=0; cur=0
for line in open(f"{d}.ont.fa"):
    if line.startswith(">"):
        if cur: tot+=cur; n+=1
        cur=0
    else: cur+=len(line.strip())
if cur: tot+=cur; n+=1
span=(106879812-105583731)+(90235370-88857161)+(22923034-22026076)
print(f"  {d}: {n} reads, {tot/1e6:.1f} Mb over {span/1e6:.2f} Mb of Ig loci -> {tot/span:.1f}x coverage")
print(f"     a clone at 10% of cells gives ~{0.5*0.10*tot/span:.1f} junction-spanning reads")
PY
done
