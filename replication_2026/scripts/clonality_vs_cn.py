#!/usr/bin/env python3
"""Correlate PoN-cleaned WGS BCR clonality (valid V+J, private in <=2 of 2,502 donors, no count<=20 cap)
with EBV CN and IR1 repeat count across full clean 1kGP (n = 2,313)."""
import csv, glob, math, os, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
cn, wt, ir1 = {}, {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try: ir1[f[0]] = float(f[5])
        except ValueError: pass

unrel = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t')}
clean = {s for s in unrel if s in cn and wt.get(s, 0) <= 0.01}

wgs_raw = {}
cdr3_donors = defaultdict(set)
for f in glob.glob(f'{BASE}/bcr/wgs/*_report.tsv'):
    s = os.path.basename(f).replace('_report.tsv', '')
    clones = []
    for line in open(f):
        if line.startswith('#') or line.startswith('count'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 8 or not p[4].startswith('IGH'): continue
        cnt = int(p[0]); cdr3 = p[3]; v = p[4]; j = p[6]
        if v == '*' or j == '*' or len(cdr3) < 8 or not cdr3.startswith('C'): continue
        clones.append((cnt, cdr3))
        cdr3_donors[cdr3].add(s)
    wgs_raw[s] = clones

public_pon = {c for c, ds in cdr3_donors.items() if len(ds) > 2}
top_frac = {}
for s, clones in wgs_raw.items():
    if s not in clean: continue
    priv = [cnt for cnt, c in clones if c not in public_pon]
    if priv and sum(priv) > 0:
        top_frac[s] = max(priv) / sum(priv)

def spearman(x, y):
    n = len(x)
    rx = [0]*n; ry = [0]*n
    for idx, (_, i) in enumerate(sorted((v, i) for i, v in enumerate(x))): rx[i] = idx
    for idx, (_, i) in enumerate(sorted((v, i) for i, v in enumerate(y))): ry[i] = idx
    mx, my = (n-1)/2.0, (n-1)/2.0
    return sum((a-mx)*(b-my) for a,b in zip(rx,ry)) / math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))

sc = [s for s in top_frac if s in cn]
si = [s for s in top_frac if s in ir1]
print(f'Clean 1kGP lines with PoN-cleaned WGS IGH clonotypes: n = {len(sc)}, median top_clone_frac = {st.median([top_frac[s] for s in sc]):.4f}')
print(f'  Spearman rho(top_clone_frac, log2_CN) = {spearman([top_frac[s] for s in sc], [math.log2(cn[s]) for s in sc]):+.4f} (n={len(sc)})')
print(f'  Spearman rho(top_clone_frac, IR1)     = {spearman([top_frac[s] for s in si], [ir1[s] for s in si]):+.4f} (n={len(si)})')
