"""Is the Europe-vs-rest IR1 split technical?

Tests the IR1 unit estimate against library-level metrics that would drive a coverage
artifact: GC bias (measured earlier on 165 samples), autosomal depth, and EBV copy number.
A real two-stock difference should be independent of all three.
"""
import csv, glob, math, os, statistics as st
from collections import defaultdict
BASE='/work/nfp8/EBV'; REF_UNITS=22737/3072.0

ir1={}
for f in glob.glob(f'{BASE}/ir1/*.tsv'):
    r=list(csv.DictReader(open(f), delimiter='\t'))
    if r:
        try:
            v=float(r[0]['ir1_copies_vs_unique'])
            if v>0: ir1[r[0]['sample']]=v*REF_UNITS
        except (ValueError,KeyError): pass
meta={}
for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t'):
    try: meta[r['sample']]=(float(r['wt_fraction']), float(r['ebv_cn_depth_ratio']),
                            float(r['auto_depth']), float(r['ebv_depth']))
    except ValueError: pass
pm={m['sample_kgpID']:(m['population'],m['continentalGroup']) for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}

gc={}
with open(f'{BASE}/metadata/window_gc.tsv') as fh:
    for line in fh:
        reg,v=line.rsplit('\t',1); c,rng=reg.split(':'); s=int(rng.split('-')[0])-1
        gc[f'{c}_{s//100000}']=float(v)
lo=[k for k,_ in sorted(gc.items(), key=lambda kv: kv[1])[:30]]
hi=[k for k,_ in sorted(gc.items(), key=lambda kv: kv[1])[-30:]]
gcb={}
for f in glob.glob(f'{BASE}/denom/*.windows.tsv'):
    s=os.path.basename(f).replace('.windows.tsv','')
    d={}
    for line in open(f):
        k,n,dep=line.split('\t')
        if int(n)>50000: d[k]=float(dep)
    h=[d[k] for k in hi if k in d]; l=[d[k] for k in lo if k in d]
    if h and l: gcb[s]=st.mean(h)/st.mean(l)

clean={s:v for s,v in ir1.items() if s in meta and meta[s][0]<=0.01}
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
    den=math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))
    return sum(p*q for p,q in zip(da,db))/den if den else 0.0
def pv(r,n):
    if n<5 or abs(r)>=1: return 0.0
    t=abs(r)*math.sqrt((n-2)/(1-r*r)); return math.erfc(t/math.sqrt(2))

print(f'\nIR1 unit estimate vs library metrics (n = {len(clean)} clean lines):')
for name, getter, subset in (
    ('autosomal depth',   lambda s: meta[s][2], clean),
    ('EBV copy number',   lambda s: math.log2(meta[s][1]), clean),
    ('EBV chrEBV depth',  lambda s: math.log2(meta[s][3]), clean),
    ('GC-bias metric',    lambda s: gcb[s], {s:v for s,v in clean.items() if s in gcb}),
):
    ss=sorted(subset); 
    if len(ss)<20: print(f'  {name:<20} n={len(ss)} too few'); continue
    r=sp([clean[s] for s in ss],[getter(s) for s in ss])
    print(f'  {name:<20} rho = {r:+.3f}  p = {pv(r,len(ss)):.1e}  (n={len(ss)})')

# Is the EUR/non-EUR split still there after regressing out autosomal depth?
ss=[s for s in sorted(clean) if s in pm]
y=[clean[s] for s in ss]; x=[meta[s][2] for s in ss]
mx,my=st.mean(x),st.mean(y)
b=sum((a-mx)*(c-my) for a,c in zip(x,y))/sum((a-mx)**2 for a in x)
res=[c-(my+b*(a-mx)) for a,c in zip(x,y)]
eur=[r for s,r in zip(ss,res) if pm[s][1]=='EUR']
oth=[r for s,r in zip(ss,res) if pm[s][1]!='EUR']
pooled=math.sqrt((st.pvariance(eur)*len(eur)+st.pvariance(oth)*len(oth))/(len(eur)+len(oth)))
print(f'\nEUR vs non-EUR, after regressing out autosomal depth:')
print(f'  EUR n={len(eur)} median residual {st.median(eur):+.3f}u;  non-EUR n={len(oth)} {st.median(oth):+.3f}u')
print(f'  Cohen d = {(st.mean(eur)-st.mean(oth))/pooled:+.2f}')
