#!/usr/bin/env bash
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
$SAM depth -a -r chrEBV:137000-154000 "$1" | awk '
{b=int($2/500); s[b]+=$3; n[b]++}
END{for(i in s) printf "%d\t%.1f\n", i*500, s[i]/n[i]}' | sort -n
