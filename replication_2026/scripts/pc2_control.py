"""Control: is the EBV-CN ~ genotype-PC2 association real, or denominator bias?

CN = 2 * ebv_depth / auto_depth, and auto_depth comes from 24 fixed GRCh38 windows.
If sequencing depth or reference bias varies with ancestry, auto_depth alone would
correlate with PC2 and manufacture the signal. Test each term separately.
"""
import csv, glob, gzip, math, statistics as st

vals = {}
for f in glob.glob('/work/nfp8/EBV/cn/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if not r:
        continue
    d = r[0]
    try:
        vals[d['sample']] = (float(d['ebv_depth']), float(d['auto_depth']),
                             float(d['ebv_cn_depth_ratio']), float(d['ebv_cn_read_ratio']))
    except ValueError:
        pass

Z = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data'
with gzip.open(f'{Z}/QTL_results/eQTL_results/expression_quants/eQTL_covariates.tab.gz', 'rt') as fh:
    rows = [l.rstrip('\n').split('\t') for l in fh]
ids, cov = rows[0][1:], {r[0]: r[1:] for r in rows[1:]}
idx = [i for i, s in enumerate(ids) if s in vals]

def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
    while i < len(o):
        j = i
        while j+1 < len(o) and v[o[j+1]] == v[o[i]]: j += 1
        for k in range(i, j+1): r[o[k]] = (i+j)/2+1
        i = j+1
    return r

def sp(a, b):
    a, b = rank(a), rank(b)
    ma, mb = st.mean(a), st.mean(b)
    da, db = [p-ma for p in a], [q-mb for q in b]
    den = math.sqrt(sum(p*p for p in da)*sum(q*q for q in db))
    return sum(p*q for p, q in zip(da, db))/den if den else 0.0

def p_of(r, n):
    t = abs(r)*math.sqrt((n-2)/(1-r*r)) if abs(r) < 1 else 99
    return math.erfc(t/math.sqrt(2))

n = len(idx)
terms = {
    'log2 EBV CN (depth ratio)': [math.log2(vals[ids[i]][2]) for i in idx],
    'log2 EBV CN (read ratio)':  [math.log2(vals[ids[i]][3]) for i in idx],
    'log2 ebv_depth (numerator)':[math.log2(vals[ids[i]][0]) for i in idx],
    'auto_depth (DENOMINATOR)':  [vals[ids[i]][1] for i in idx],
}
print(f'n = {n} MAGE donors\n')
print(f'{"term":<28} {"vs PC1":>9} {"vs PC2":>9} {"p(PC2)":>10} {"vs PC3":>9}')
for name, v in terms.items():
    rs = {}
    for pc in ('PC1', 'PC2', 'PC3'):
        y = [float(cov[pc][i]) for i in idx]
        rs[pc] = sp(v, y)
    print(f'{name:<28} {rs["PC1"]:>9.3f} {rs["PC2"]:>9.3f} {p_of(rs["PC2"], n):>10.2e} {rs["PC3"]:>9.3f}')

# the two estimators share a denominator; do they agree on the PC2 signal?
print('\nestimator concordance (Spearman between the two CN estimates): '
      f'{sp(terms["log2 EBV CN (depth ratio)"], terms["log2 EBV CN (read ratio)"]):.4f}')
