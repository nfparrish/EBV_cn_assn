"""Decisive control: does GC bias explain the CN population structure?

GC bias is itself population-structured (eta-sq 0.091). Regress log2 CN on the GC-bias
metric and re-test the population structure of the residual. If the signal survives,
GC bias is not the driver.
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

ebvd = {}
with open('/work/nfp8/EBV/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try: ebvd[r['sample']] = float(r['ebv_depth'])
        except ValueError: pass
pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}

S, Y, G = [], [], []
for f in glob.glob('/work/nfp8/EBV/denom/*.windows.tsv'):
    s = os.path.basename(f).replace('.windows.tsv', '')
    if s not in ebvd or s not in pop: continue
    d = {}
    with open(f) as fh:
        for line in fh:
            k, n, dep = line.split('\t')
            if int(n) > 50000: d[k] = float(dep)
    if len(d) < 140: continue
    hi = [d[k] for k in hi_k if k in d]; lo = [d[k] for k in lo_k if k in d]
    if not hi or not lo: continue
    S.append(s); Y.append(math.log2(2*ebvd[s]/st.median(d.values()))); G.append(st.mean(hi)/st.mean(lo))

n = len(S)
mg, my = st.mean(G), st.mean(Y)
b = sum((g-mg)*(y-my) for g, y in zip(G, Y)) / sum((g-mg)**2 for g in G)
R = [y - (my + b*(g-mg)) for y, g in zip(Y, G)]
print(f'n = {n};  slope of log2 CN on GC-bias = {b:+.2f}  '
      f'(variance explained {100*b*b*st.pvariance(G)/st.pvariance(Y):.1f}%)\n')

def rank(v):
    o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
    while i<len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r

def eta(vals):
    g = defaultdict(list)
    for s, v in zip(S, vals): g[pop[s]].append(v)
    g = {k: v for k, v in g.items() if len(v) >= 5}
    flat=[v for gg in g.values() for v in gg]; rk=rank(flat); pos=0; H=0.0
    for gg in g.values():
        rs=sum(rk[pos:pos+len(gg)]); pos+=len(gg); H += rs*rs/len(gg)
    H = 12.0/(len(flat)*(len(flat)+1))*H - 3*(len(flat)+1); k=len(g)
    return max(0.0,(H-k+1)/(len(flat)-k)), k, len(flat)

e1, k, nn = eta(Y); e2, _, _ = eta(R); e3, _, _ = eta(G)
print(f'population structure (eta-sq on ranks, {k} populations, n={nn}):')
print(f'  log2 CN, raw                  = {e1:.4f}')
print(f'  log2 CN, GC-bias residualised = {e2:.4f}   <-- survives? {"YES" if e2 > 0.6*e1 else "NO"}')
print(f'  GC-bias metric itself         = {e3:.4f}')
print(f'\nretained after control: {100*e2/e1:.0f}% of the raw population signal')
