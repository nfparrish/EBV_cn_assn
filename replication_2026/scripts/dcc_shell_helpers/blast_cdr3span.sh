#!/usr/bin/env bash
# Does requiring the alignment to SPAN THE CDR3 restore clone specificity?
# Queries are CDR3 -120/+80, so the CDR3 occupies query positions 121..121+cdr3len.
# A clone-specific hit must cover that interval; a germline-V-driven hit will not.
set -euo pipefail
cd /work/nfp8/EBV/bcr_blast/proto
BL=/hpc/group/parrishlab/envs/vcall-env/bin
for d in HG00096 HG00130; do
  tag=$([ "$d" = HG00096 ] && echo OWN || echo CTRL)
  $BL/blastn -query HG00096_junction_queries.fa -db ${d}_db \
     -outfmt '6 qseqid sseqid pident length qlen qstart qend' \
     -evalue 1e-5 -word_size 11 -num_threads 2 2>/dev/null \
  | awk -v TAG="$tag" -v D="$d" '
      { split($1,f,"|"); cl=0
        for (i in f) if (f[i] ~ /^cdr3len=/) { split(f[i],g,"="); cl=g[2] }
        cs=121; ce=120+cl                      # CDR3 interval within the query
        anyhit[$1]=1
        if ($3>=85 && $6<=cs+5 && $7>=ce-5) { span[$1]=1; nspan++ } }
      END { nq=0; for (q in anyhit) nq++
            ns=0; for (q in span) ns++
            printf "  %-6s(%-4s): queries with ANY hit %2d   queries with CDR3-SPANNING hit %2d   spanning alignments %3d\n",
                   D, TAG, nq, ns, nspan+0 }'
done
