#!/usr/bin/env bash
#SBATCH --job-name=orilyt_chk
#SBATCH --account=parrishlab
#SBATCH --partition=common,scavenger
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/work/nfp8/EBV/logs/orilyt_chk.out
/hpc/group/parrishlab/envs/pysam/bin/python3 -u - << 'PYEOF'
import csv, os, pysam, statistics as st
from concurrent.futures import ProcessPoolExecutor

BASE = '/work/nfp8/EBV'
cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try:
            v = float(f[3]); w = float(f[9])
            if v > 0: cn[f[0]] = v; wt[f[0]] = w
        except ValueError: pass

unrel = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t') if r['sample'] in cn}
clean = sorted([s for s in unrel if wt[s] <= 0.01], key=lambda s: cn[s])
low30  = [s for s in clean[:35] if os.path.exists(f'{BASE}/slices/{s}.chrEBV.bam')][:30]
high30 = [s for s in clean[-35:] if os.path.exists(f'{BASE}/slices/{s}.chrEBV.bam')][:30]

def get_metrics(s):
    bam = pysam.AlignmentFile(f'{BASE}/slices/{s}.chrEBV.bam', 'rb')
    c_ori = bam.count('chrEBV', 41000, 55000) / 14000.0
    c_ant = bam.count('chrEBV', 110000, 130000) / 20000.0
    tr_tot = tr_clip = 0
    for r in bam.fetch('chrEBV', 169000, 171823):
        if r.is_unmapped or r.is_secondary: continue
        tr_tot += 1
        if r.cigartuples and ((r.cigartuples[0][0] == 4 and r.cigartuples[0][1] >= 15) or (r.cigartuples[-1][0] == 4 and r.cigartuples[-1][1] >= 15)):
            tr_clip += 1
    bam.close()
    return (c_ori / c_ant if c_ant > 0 else 1.0), (tr_clip / tr_tot if tr_tot > 0 else 0.0)

with ProcessPoolExecutor(max_workers=4) as ex:
    m_low  = list(ex.map(get_metrics, low30))
    m_high = list(ex.map(get_metrics, high30))

print(f'=== 9. oriLyt (41-55kb) / Antipode (110-130kb) & TR CLEAVAGE ===')
print(f'  Bottom 30 low-CN  (median CN {st.median([cn[s] for s in low30]):5.1f}): oriLyt/antipode median = {st.median([x[0] for x in m_low]):.4f} (mean {st.mean([x[0] for x in m_low]):.4f}), TR clip frac = {st.median([x[1] for x in m_low]):.4f}')
print(f'  Top 30 high-CN    (median CN {st.median([cn[s] for s in high30]):5.1f}): oriLyt/antipode median = {st.median([x[0] for x in m_high]):.4f} (mean {st.mean([x[0] for x in m_high]):.4f}), TR clip frac = {st.median([x[1] for x in m_high]):.4f}')
PYEOF
