"""Three-platform EBV CN comparison, and a check that downsampling defused the
CN-dependent variant-detection confound.
"""
import csv, glob, math, os, statistics as st

def load_alt(tag, field='ebv_cn'):
    d = {}
    for f in glob.glob(f'/work/nfp8/EBV/{tag}/*.tsv'):
        with open(f) as fh:
            r = list(csv.DictReader(fh, delimiter='\t'))
        if r:
            try:
                v = float(r[0][field]); ed = float(r[0]['ebv_depth'])
                if v > 0 and ed > 2.0: d[r[0]['sample']] = v
            except (ValueError, KeyError): pass
    return d

short = {}
with open('/work/nfp8/EBV/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            v = float(r['ebv_cn_depth_ratio'])
            if v > 0: short[r['sample']] = v
        except ValueError: pass
ont, low = load_alt('ont'), load_alt('lowcov')

def rank(v):
    o=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
    while i<len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r
def pear(a,b):
    ma,mb=st.mean(a),st.mean(b); da=[p-ma for p in a]; db=[q-mb for q in b]
    return sum(p*q for p,q in zip(da,db))/math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))

print('=== platform comparison (log2 CN) ===')
for name, d in (('ONT (~2023, long-read)', ont), ('low-coverage (2008-12)', low)):
    sh = sorted(set(d) & set(short))
    if len(sh) < 20: print(f'{name}: only {len(sh)} overlapping'); continue
    x = [math.log2(short[s]) for s in sh]; y = [math.log2(d[s]) for s in sh]
    ratios = [d[s]/short[s] for s in sh]; lr = [math.log2(v) for v in ratios]
    print(f'\n{name}  n={len(sh)}')
    print(f'  vs 30x: Pearson r={pear(x,y):.3f}  Spearman={pear(rank(x),rank(y)):.3f}  R2={pear(x,y)**2:.3f}')
    print(f'  ratio to 30x: median {st.median(ratios):.3f}  geo-mean {2**st.mean(lr):.3f}  '
          f'sd(log2) {st.pstdev(lr):.3f}')

both = sorted(set(ont) & set(low) & set(short))
if len(both) >= 20:
    print(f'\n=== all three platforms, n={len(both)} ===')
    o = [math.log2(ont[s]) for s in both]; l = [math.log2(low[s]) for s in both]
    s3 = [math.log2(short[s]) for s in both]
    print(f'  30x vs ONT      r={pear(s3,o):.3f}')
    print(f'  30x vs lowcov   r={pear(s3,l):.3f}')
    print(f'  ONT vs lowcov   r={pear(o,l):.3f}   <-- two independent non-NYGC measurements')

print('\n=== did downsampling defuse the CN-detection confound? ===')
rows = []
with open('/work/nfp8/EBV/variants/lofreq_summary.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4:
            try: rows.append((p[0], float(p[1]), float(p[2]), int(p[3])))
            except ValueError: pass
print(f'n = {len(rows)} samples called')
raw = [math.log2(r[1]) for r in rows]; eff = [r[2] for r in rows]; nv = [r[3] for r in rows]
print(f'  effective depth: mean {st.mean(eff):.1f}x  sd {st.pstdev(eff):.2f}  '
      f'range {min(eff):.1f}-{max(eff):.1f}')
print(f'  variants/sample: mean {st.mean(nv):.1f}  median {st.median(nv):.0f}  range {min(nv)}-{max(nv)}')
print(f'  n_variants vs log2(RAW depth): r = {pear(raw, nv):+.3f}   <-- want ~0')
cn = [math.log2(short[r[0]]) for r in rows if r[0] in short]
nv2 = [r[3] for r in rows if r[0] in short]
print(f'  n_variants vs log2(CN)       : r = {pear(cn, nv2):+.3f}   <-- want ~0 (n={len(cn)})')
