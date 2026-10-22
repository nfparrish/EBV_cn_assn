#!/usr/bin/env bash
SAM=/hpc/group/parrishlab/envs/samtools-env/bin/samtools
$SAM depth -a -r chrEBV "$1" | awk '
{b=int($2/4000); s[b]+=$3; n[b]++}
END{for(i=0;i<=42;i++){d=(n[i]?s[i]/n[i]:0); printf "%6d-%6d %9.1f %s\n", i*4000, (i+1)*4000, d, (d<20?"<== LOW":"")}}'
