"""Does the per-sample GC-bias metric confound EBV CN?

EBV is 59.4% GC, the autosomal panel 40.8%. If GC bias varied across samples it would
act differentially on numerator and denominator and could manufacture CN differences --
including population-structured ones. Test directly.
"""
import csv, glob, math, os, statistics as st
from collections import defaultdict

gc = {}
with open('/work/nfp8/EBV/metadata/window_gc.tsv') as fh:
    for line in fh:
        reg, v = line.rsplit('\t', 1)
        c, rng = reg.split(':'); s = int(rng.split('-')[0]) - 1
        gc[f'{c}_{s//100000}'] = float(v)
lo_k = [k for k, _ in sorted(gc.items(), key=lambda kv: kv[1])[:30]]
hi_k = [k for k, _ in sorted(gc.items(), key=lambda kv: kv[1])[-30:]]

old = {}
with open('/work/nfp8/EBV/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try: old[r['sample']] = (float(r['ebv_cn_depth_ratio']), float(r['ebv_depth']))
        except ValueError: pass

rows = []
for f in glob.glob('/work/nfp8/EBV/denom/*.windows.tsv'):
    s = os.path.basename(f).replace('.windows.tsv', '')
    if s not in old: continue
    d = {}
    with open(f) as fh:
        for line in fh:
            k, n, dep = line.split('\t')
            if int(n) > 50000: d[k] = float(dep)
    if len(d) < 140: continue
    hi = [d[k] for k in hi_k if k in d]; lo = [d[k] for k in lo_k if k in d]
    if not hi or not lo: continue
    rows.append((s, st.median(d.values()), st.mean(hi)/st.mean(lo), old[s][0], old[s][1]))

print(f'n = {len(rows)} samples\n')
def rank(v):
    o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
    while i<len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r
def sp(a,b):
    a,b=rank(a),rank(b); ma,mb=st.mean(a),st.mean(b)
    da=[p-ma for p in a]; db=[q-mb for q in b]
    return sum(p*q for p,q in zip(da,db))/math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))
def pv(r,n):
    t=abs(r)*math.sqrt((n-2)/(1-r*r)) if abs(r)<1 else 99
    return math.erfc(t/math.sqrt(2))

gcb = [r[2] for r in rows]; cn_new = [math.log2(2*r[4]/r[1]) for r in rows]
cn_old = [math.log2(r[3]) for r in rows]; n = len(rows)
print(f'GC-bias metric: mean {st.mean(gcb):.3f}  sd {st.pstdev(gcb):.3f}  range {min(gcb):.3f}-{max(gcb):.3f}')
print(f'\n  GC-bias vs log2 CN        : rho = {sp(gcb, cn_new):+.3f}  p = {pv(sp(gcb,cn_new),n):.3f}')
print(f'  GC-bias vs autosomal depth: rho = {sp(gcb, [r[1] for r in rows]):+.3f}')
print(f'  CN_old vs CN_new          : rho = {sp(cn_old, cn_new):.5f}')

pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}
g = defaultdict(list)
for r in rows:
    if r[0] in pop: g[pop[r[0]]].append(r[2])
g = {k: v for k, v in g.items() if len(v) >= 5}
if len(g) > 2:
    flat=[v for gg in g.values() for v in gg]; rk=rank(flat); pos=0; H=0.0
    for gg in g.values():
        rs=sum(rk[pos:pos+len(gg)]); pos+=len(gg); H += rs*rs/len(gg)
    H = 12.0/(len(flat)*(len(flat)+1))*H - 3*(len(flat)+1); k=len(g)
    print(f'\n  GC-bias population structure: eta-sq = {max(0,(H-k+1)/(len(flat)-k)):.4f} over {k} populations')
    print('  (compare: CN eta-sq = 0.126, denominator eta-sq = 0.017)')
