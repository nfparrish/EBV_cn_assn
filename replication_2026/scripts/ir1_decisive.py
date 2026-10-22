#!/usr/bin/env python3
"""Full-cohort IR1 repeat validation: within-chrEBV GC bias control (n=2,449) and ONT long-read cross-platform check (n=891)."""
import csv, math, statistics as st

BASE = '/work/nfp8/EBV'
ped = {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    if line.startswith('FamilyID'): continue
    f = line.rstrip('\n').split('\t')
    ped[f[1]] = (f[5], f[6]) # pop, superpop

wt = {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: wt[f[0]] = float(f[9])
        except ValueError: pass

ir1 = {}
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try: ir1[f[0]] = float(f[5])
        except ValueError: pass

gcb = {}
for line in open(f'{BASE}/metadata/all_gcbias_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5 and f[4] != 'NA':
        try: gcb[f[0]] = float(f[4])
        except ValueError: pass

ont_ir1 = {}
for line in open(f'{BASE}/metadata/all_ont_ir1_907.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5:
        try: ont_ir1[f[0]] = float(f[4])
        except ValueError: pass

unrel = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t')}
clean = [s for s in unrel if s in ir1 and wt.get(s, 0) <= 0.01 and s in ped]

def cohens_d(a, b):
    ma, mb = st.mean(a), st.mean(b)
    sp = math.sqrt(((len(a)-1)*st.variance(a) + (len(b)-1)*st.variance(b)) / (len(a)+len(b)-2))
    return (ma - mb) / sp

# 1. GC bias residualization at full n = 2,449
sgc = [s for s in clean if s in gcb]
x = [gcb[s] for s in sgc]; y = [ir1[s] for s in sgc]
mx, my = st.mean(x), st.mean(y)
b = sum((xi-mx)*(yi-my) for xi,yi in zip(x,y)) / sum((xi-mx)**2 for xi in x)
res = {s: ir1[s] - (my + b*(gcb[s]-mx)) for s in sgc}
eur4 = {'CEU', 'FIN', 'IBS', 'GBR'}
d_raw = cohens_d([ir1[s] for s in sgc if ped[s][0] in eur4], [ir1[s] for s in sgc if ped[s][0] not in eur4])
d_adj = cohens_d([res[s] for s in sgc if ped[s][0] in eur4], [res[s] for s in sgc if ped[s][0] not in eur4])
print(f'1. Within-chrEBV GC control (n={len(sgc)}): raw d = {d_raw:.3f} -> GC-adjusted d = {d_adj:.3f} ({d_adj/d_raw*100:.1f}% survives)')

# 2. ONT cross-platform validation at full n = 891
sont = [s for s in clean if s in ont_ir1]
xo = [ir1[s] for s in sont]; yo = [ont_ir1[s] for s in sont]
mxo, myo = st.mean(xo), st.mean(yo)
r = sum((a-mxo)*(b-myo) for a,b in zip(xo,yo)) / math.sqrt(sum((a-mxo)**2 for a in xo)*sum((b-myo)**2 for b in yo))
d_ont = cohens_d([ont_ir1[s] for s in sont if ped[s][0] in eur4], [ont_ir1[s] for s in sont if ped[s][0] not in eur4])
print(f'2. ONT vs 30x IR1 (n={len(sont)}): r = {r:.4f}, median ratio ONT/30x = {st.median([yo[i]/xo[i] for i in range(len(sont))]):.4f}, ONT 4-EUR d = {d_ont:.3f}')
