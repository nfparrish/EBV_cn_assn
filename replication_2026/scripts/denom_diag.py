"""Diagnose the autosomal denominator.

The 1kGP high-coverage project targeted 30x, so auto_depth should cluster near 30 with
modest spread. Critically, the headline population result is only safe if the DENOMINATOR
carries no population structure -- otherwise CN differences could be denominator artifacts.
"""
import csv, glob, math, statistics as st
from collections import defaultdict

rows = {}
for f in glob.glob('/work/nfp8/EBV/cn/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            rows[r[0]['sample']] = {k: float(r[0][k]) for k in
                ('ebv_depth','auto_depth','ebv_cn_depth_ratio','auto_reads','read_len')}
        except ValueError:
            pass
print(f'n = {len(rows)} samples\n')

ad = sorted(v['auto_depth'] for v in rows.values())
def q(p): return ad[int(len(ad)*p)]
print('autosomal depth (expected ~30x for this project):')
print(f'  min {ad[0]:.1f}  p1 {q(.01):.1f}  p25 {q(.25):.1f}  median {q(.5):.1f} '
      f' p75 {q(.75):.1f}  p99 {q(.99):.1f}  max {ad[-1]:.1f}')
print(f'  mean {st.mean(ad):.1f}  sd {st.pstdev(ad):.1f}  CV {st.pstdev(ad)/st.mean(ad):.3f}')
lo = [s for s, v in rows.items() if v['auto_depth'] < 20]
hi = [s for s, v in rows.items() if v['auto_depth'] > 60]
print(f'  <20x: {len(lo)}   >60x: {len(hi)}')

rl = set(round(v['read_len']) for v in rows.values())
print(f'\nread lengths observed: {sorted(rl)}')

# population structure of the DENOMINATOR -- the key control
pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}

def rank_all(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r=[0.0]*len(v); i=0
    while i < len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r

def kruskal_eta(groups):
    flat=[v for g in groups.values() for v in g]; n=len(flat)
    rk=rank_all(flat); pos=0; H=0.0
    for g in groups.values():
        rs=sum(rk[pos:pos+len(g)]); pos+=len(g); H += rs*rs/len(g)
    H = 12.0/(n*(n+1))*H - 3*(n+1); k=len(groups)
    return H, k, max(0.0,(H-k+1)/(n-k))

print('\npopulation structure, by term (eta-sq on ranks, 26 populations):')
for term, label in (('auto_depth','DENOMINATOR auto_depth'),
                    ('ebv_depth','numerator ebv_depth'),
                    ('ebv_cn_depth_ratio','CN (ratio)')):
    g = defaultdict(list)
    for s, v in rows.items():
        if s in pop:
            g[pop[s]].append(math.log2(v[term]) if term != 'auto_depth' else v[term])
    g = {k2: v for k2, v in g.items() if len(v) >= 20}
    H, k, eta = kruskal_eta(g)
    print(f'  {label:<24} eta-sq = {eta:.4f}   H = {H:7.1f}')
