#!/usr/bin/env python3
"""Population and continental structure of IR1 repeat copy number at full n = 2,449 clean 1kGP."""
import csv, math, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
ped = {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    if line.startswith('FamilyID'): continue
    f = line.rstrip('\n').split('\t')
    ped[f[1]] = (f[5], f[6]) # pop, superpop

cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass

ir1 = {}
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try: ir1[f[0]] = float(f[5])
        except ValueError: pass

unrel = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t')}
clean = [s for s in unrel if s in ir1 and wt.get(s, 0) <= 0.01 and s in ped]
print(f'Clean unrelated 1kGP with IR1: n = {len(clean)}, median = {st.median([ir1[s] for s in clean]):.2f}u')
by_pop = defaultdict(list)
for s in clean:
    by_pop[ped[s][0]].append(ir1[s])
for p in sorted(by_pop, key=lambda k: st.median(by_pop[k])):
    print(f'  {p:4s} ({ped[[s for s in clean if ped[s][0]==p][0]][1]}): n={len(by_pop[p]):3d}, median={st.median(by_pop[p]):5.2f}u, mean={st.mean(by_pop[p]):5.2f}u')
