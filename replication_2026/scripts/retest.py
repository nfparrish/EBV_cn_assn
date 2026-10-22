"""Test-retest of EBV CN: NYGC 30x (2019) vs phase 3 low-coverage (2008-2012 reads).

Same individuals, different passage cells, ~decade apart, independent libraries and
sequencing centres. Reproducibility here separates "stable line property" from
"culture noise" -- though NOT from "set at transformation then clonally fixed",
which looks identical and carries zero heritability.
"""
import csv, glob, math, os, statistics as st

new = {}
for f in glob.glob('/work/nfp8/EBV/lowcov/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            v = float(r[0]['ebv_cn']); d = float(r[0]['auto_median_depth'])
            if v > 0: new[r[0]['sample']] = (v, d)
        except (ValueError, KeyError): pass

old = {}
with open('/work/nfp8/EBV/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            v = float(r['ebv_cn_depth_ratio'])
            if v > 0: old[r['sample']] = v
        except (ValueError, KeyError): pass

both = sorted(set(new) & set(old))
print(f'n = {len(both)} samples with both measurements\n')
x = [math.log2(old[s]) for s in both]
y = [math.log2(new[s][0]) for s in both]

def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r=[0.0]*len(v); i=0
    while i < len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r

def pear(a,b):
    ma,mb=st.mean(a),st.mean(b); da=[p-ma for p in a]; db=[q-mb for q in b]
    return sum(p*q for p,q in zip(da,db))/math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))

r_p, r_s = pear(x,y), pear(rank(x),rank(y))
ratios = [new[s][0]/old[s] for s in both]
lr = [math.log2(v) for v in ratios]
print(f'Pearson r (log2 CN)   = {r_p:.4f}')
print(f'Spearman rho          = {r_s:.4f}')
print(f'R-squared             = {r_p**2:.4f}')
print(f'\nlowcov/30x ratio: median {st.median(ratios):.3f}  geometric mean {2**st.mean(lr):.3f}')
print(f'  sd of log2 ratio = {st.pstdev(lr):.3f}  (= {100*(2**st.pstdev(lr)-1):.1f}% typical deviation)')
within = sum(1 for v in ratios if 0.5 <= v <= 2.0)
print(f'  within 2-fold: {within}/{len(ratios)} ({100*within/len(ratios):.0f}%)')
ad = [new[s][1] for s in both]
print(f'\nlow-coverage autosomal depth: median {st.median(ad):.1f}x  range {min(ad):.1f}-{max(ad):.1f}x')
print(f'  (depth does NOT correlate with the ratio: r = {pear([math.log2(v) for v in ad], lr):+.3f})')
print('\nlargest disagreements:')
for s in sorted(both, key=lambda s: -abs(math.log2(new[s][0]/old[s])))[:6]:
    print(f'  {s:<10} 30x={old[s]:>8.2f}  lowcov={new[s][0]:>8.2f}  ratio={new[s][0]/old[s]:.2f}  (lowcov depth {new[s][1]:.1f}x)')
