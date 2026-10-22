"""Does EBV CN structure with continental group, population, or MAGE sequencing batch?

Batch is a shipment-and-sequencing batch, so it is the only culture/handling proxy in
the metadata. Ancestry and batch may themselves be confounded (1kGP populations were
collected and transformed at different times), which is the crux of the
individual-vs-culture aliasing problem.
Kruskal-Wallis on log2 CN; variance explained via eta-squared on ranks.
"""
import csv, glob, math, statistics as st
from collections import defaultdict

cn = {}
for f in glob.glob('/work/nfp8/EBV/cn/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            v = float(r[0]['ebv_cn_depth_ratio'])
            if v > 0: cn[r[0]['sample']] = math.log2(v)
        except ValueError: pass

meta = list(csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t'))
rows = [(m, cn[m['sample_kgpID']]) for m in meta if m['sample_kgpID'] in cn]
print(f'{len(rows)} MAGE QTL-set donors with EBV CN\n')

def rank_all(vals):
    o = sorted(range(len(vals)), key=lambda i: vals[i]); r=[0.0]*len(vals); i=0
    while i < len(o):
        j=i
        while j+1<len(o) and vals[o[j+1]]==vals[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r

def kruskal(groups):
    flat=[v for g in groups.values() for v in g]; n=len(flat)
    rk=rank_all(flat); pos=0; H=0.0
    for g in groups.values():
        rs=sum(rk[pos:pos+len(g)]); pos+=len(g)
        H += rs*rs/len(g)
    H = 12.0/(n*(n+1))*H - 3*(n+1)
    k=len(groups)
    # eta^2 on ranks, and a chi2 tail via Wilson-Hilferty
    eta = max(0.0, (H-k+1)/(n-k))
    x=H/(k-1); p=math.erfc(((x**(1/3)-(1-2/(9*(k-1))))*math.sqrt(9*(k-1)/2))/math.sqrt(2))/1 if k>1 else 1
    return H, k, eta, min(1.0,p)

for field, label in (('continentalGroup','continental group'), ('population','population'), ('batch','batch')):
    g = defaultdict(list)
    for m, v in rows: g[m[field]].append(v)
    g = {k: v for k, v in g.items() if len(v) >= 5}
    H, k, eta, p = kruskal(g)
    print(f'{label:<20} groups={k:<3} H={H:7.1f}  eta^2(ranks)={eta:.3f}  p~{p:.1e}')
    med = sorted(((st.median(v), k2, len(v)) for k2, v in g.items()), reverse=True)
    print('   highest:', ', '.join(f'{k2} {2**m:.0f}cp (n={n})' for m, k2, n in med[:3]))
    print('   lowest: ', ', '.join(f'{k2} {2**m:.0f}cp (n={n})' for m, k2, n in med[-3:]))
    print()
