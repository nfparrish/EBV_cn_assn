"""IR1 (BamHI-W) repeat-unit count: distribution, population structure, and joint test
with the co-segregation lineages.

NC_007605 contains the FULL IR1 array, not one collapsed unit, so depth(IR1)/depth(unique)
estimates sample_units / reference_units. Reference array = 22,737 bp at 3,072 bp per W unit
= 7.40 units, so units = ratio * 7.40.

J Virol 2017 reports IR1 unit count differing BETWEEN B95-8 preparations (11.6 in the B95-8
cell line vs 6.6 in the HB9 BAC), which makes unit count a documented stock-level difference
and therefore an independent provenance marker -- in the one region our SNV pipeline masks.
"""
import csv, glob, math, os, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
REF_UNITS = 22737 / 3072.0

ir1 = {}
for f in glob.glob(f'{BASE}/ir1/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            v = float(r[0]['ir1_copies_vs_unique'])
            if v > 0: ir1[r[0]['sample']] = v * REF_UNITS
        except (ValueError, KeyError): pass

wt, cn = {}, {}
with open(f'{BASE}/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            wt[r['sample']] = float(r['wt_fraction']); cn[r['sample']] = float(r['ebv_cn_depth_ratio'])
        except ValueError: pass
pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/'
         'sample.metadata.MAGE.v1.0.txt'), delimiter='\t')}
coinf = {s for s, v in wt.items() if v > 0.01}
clean = {s: v for s, v in ir1.items() if s not in coinf}

v = sorted(clean.values())
def q(p): return v[int(len(v)*p)]
print(f'IR1 unit count, {len(v)} clean lines (reference array = {REF_UNITS:.2f} units)')
print(f'  min {v[0]:.1f}  p10 {q(.10):.1f}  p25 {q(.25):.1f}  median {q(.50):.1f} '
      f' p75 {q(.75):.1f}  p90 {q(.90):.1f}  max {v[-1]:.1f}')
print(f'  mean {st.mean(v):.2f}  sd {st.pstdev(v):.2f}  CV {st.pstdev(v)/st.mean(v):.3f}')
print(f'  published: B95-8 cell line 11.6, HB9 BAC 6.6\n')

def rank(x):
    o=sorted(range(len(x)),key=lambda i:x[i]); r=[0.0]*len(x); i=0
    while i<len(o):
        j=i
        while j+1<len(o) and x[o[j+1]]==x[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r
def sp(a,b):
    a,b=rank(a),rank(b); ma,mb=st.mean(a),st.mean(b)
    da=[p-ma for p in a]; db=[p-mb for p in b]
    return sum(p*q for p,q in zip(da,db))/math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))

both=[s for s in clean if s in cn]
print(f'IR1 units vs EBV copy number: rho = {sp([clean[s] for s in both],[math.log2(cn[s]) for s in both]):+.3f} (n={len(both)})')

g = defaultdict(list)
for s, u in clean.items():
    if s in pop: g[pop[s]].append(u)
g = {k: x for k, x in g.items() if len(x) >= 15}
flat=[x for xs in g.values() for x in xs]; rk=rank(flat); pos=0; H=0.0
for xs in g.values():
    rs=sum(rk[pos:pos+len(xs)]); pos+=len(xs); H += rs*rs/len(xs)
H = 12.0/(len(flat)*(len(flat)+1))*H - 3*(len(flat)+1); k=len(g)
print(f'\npopulation structure of IR1 unit count: eta-sq = {max(0,(H-k+1)/(len(flat)-k)):.4f} '
      f'over {k} populations (n={len(flat)})')
print('  (compare: EBV CN eta-sq = 0.126)')
ms = sorted(((st.median(x), p, len(x)) for p, x in g.items()), reverse=True)
print('  highest: ' + ', '.join(f'{p} {m:.1f}u' for m, p, n in ms[:4]))
print('  lowest : ' + ', '.join(f'{p} {m:.1f}u' for m, p, n in ms[-4:]))

# joint test: do the GBR / IBS lineage carriers differ in IR1 unit count?
carr = defaultdict(set)
with open(f'{BASE}/variants/all_calls.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 5 and p[1] in ('132965', '91192'): carr[p[1]].add(p[0])

print('\n--- joint test: IR1 unit count by lineage carriage ---')
for site, name in (('132965', 'GBR lineage (9-site)'), ('91192', 'IBS lineage')):
    c = [clean[s] for s in carr[site] if s in clean]
    nc = [u for s, u in clean.items() if s not in carr[site]]
    if len(c) < 5: continue
    pooled = math.sqrt((st.pvariance(c)*len(c) + st.pvariance(nc)*len(nc))/(len(c)+len(nc)))
    d = (st.mean(c)-st.mean(nc))/pooled if pooled else 0
    se = pooled*math.sqrt(1/len(c)+1/len(nc))
    z = (st.mean(c)-st.mean(nc))/se if se else 0
    print(f'  {name:<22} carriers n={len(c):<4} median {st.median(c):.2f}u   '
          f'non-carriers median {st.median(nc):.2f}u   Cohen d={d:+.2f}  z={z:+.1f}')
