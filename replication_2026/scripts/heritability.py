#!/usr/bin/env python3
"""Trio heritability of EBV copy number in 1kGP LCLs (fixed: uses f[5] 26-population column
and reports both pooled and within-population demeaned parent-offspring vs spouse-pair r)."""
import math, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
ped = {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    if line.startswith('FamilyID'): continue
    f = line.rstrip('\n').split('\t')
    ped[f[1]] = (f[0], f[2], f[3], int(f[4]), f[5], f[6]) # fam, dad, mom, sex, pop(26), superpop(5)

cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try:
            v = float(f[3]); w = float(f[9])
            if v > 0: cn[f[0]] = v; wt[f[0]] = w
        except ValueError: pass

clean = {s: math.log2(v) for s, v in cn.items() if wt.get(s, 0) <= 0.01}
trios = []
for s, (fam, fa, mo, sex, pop, spop) in ped.items():
    if fa != '0' and mo != '0' and s in clean and fa in clean and mo in clean:
        trios.append((s, fa, mo, pop, spop))

def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx)**2 for a in x) * sum((b - my)**2 for b in y))
    return num / den if den > 0 else 0.0

print(f'Clean trios: n = {len(trios)}')
for mode in ('pooled', 'within_pop'):
    p_vals = defaultdict(list)
    for s, fa, mo, p, _ in trios:
        p_vals[p].extend([clean[s], clean[fa], clean[mo]])
    p_mean = {p: st.mean(v) for p, v in p_vals.items()} if mode == 'within_pop' else defaultdict(float)
    val = lambda s, p: clean[s] - p_mean[p]
    r_mid = pearson([0.5*(val(fa,p)+val(mo,p)) for s,fa,mo,p,_ in trios], [val(s,p) for s,fa,mo,p,_ in trios])
    r_sp  = pearson([val(fa,p) for s,fa,mo,p,_ in trios] + [val(mo,p) for s,fa,mo,p,_ in trios], [val(s,p) for s,fa,mo,p,_ in trios]*2)
    r_pp  = pearson([val(fa,p) for s,fa,mo,p,_ in trios], [val(mo,p) for s,fa,mo,p,_ in trios])
    print(f'  [{mode:10s}] r_mid = {r_mid:+.4f}, r_single_parent = {r_sp:+.4f}, r_spouse = {r_pp:+.4f} -> h2 = 2(r_sp - r_pp) = {2*(r_sp-r_pp):+.4f}')
