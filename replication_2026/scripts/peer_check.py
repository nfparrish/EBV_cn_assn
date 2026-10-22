"""Correlate 1kGP-derived EBV copy number against MAGE's published eQTL covariates.

Go/no-go for the "viral state is an unmodelled confounder in MAGE" section: if a PEER
factor already absorbs viral state, that section weakens a lot and we want to know now.
Spearman + Pearson on log2 CN, with a Benjamini-Hochberg correction over the 67 covariates.
"""
import csv, glob, gzip, math, os, statistics as st

EBV = '/work/nfp8/EBV/cn'
ZEN = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data'

cn = {}
for f in glob.glob(f'{EBV}/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if not r:
        continue
    row = r[0]
    try:
        v = float(row['ebv_cn_depth_ratio'])
    except (ValueError, KeyError):
        continue
    if v > 0:
        cn[row['sample']] = v
print(f'EBV CN available for {len(cn)} samples')

with gzip.open(f'{ZEN}/QTL_results/eQTL_results/expression_quants/eQTL_covariates.tab.gz', 'rt') as fh:
    rows = [l.rstrip('\n').split('\t') for l in fh]
ids, cov = rows[0][1:], {r[0]: r[1:] for r in rows[1:]}
print(f'MAGE covariates: {len(cov)} covariates x {len(ids)} samples')

shared = [i for i, s in enumerate(ids) if s in cn]
print(f'overlap with MAGE QTL set: {len(shared)} donors')
x = [math.log2(cn[ids[i]]) for i in shared]

def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(order):                      # average ties
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    da, db = [p - ma for p in a], [q - mb for q in b]
    den = math.sqrt(sum(p * p for p in da) * sum(q * q for q in db))
    return sum(p * q for p, q in zip(da, db)) / den if den else 0.0

def pval(r, n):                                 # two-sided t approximation
    if n < 5 or abs(r) >= 1:
        return 0.0
    t = abs(r) * math.sqrt((n - 2) / (1 - r * r))
    # normal approximation to the t tail, adequate at n~700
    return math.erfc(t / math.sqrt(2))

xr = rank(x)
out = []
for name, vals in cov.items():
    try:
        y = [float(vals[i]) for i in shared]
    except ValueError:
        continue                                # non-numeric (e.g. sex as XX/XY)
    if len(set(y)) < 3:
        continue
    out.append((name, pearson(x, y), pearson(xr, rank(y)), pval(pearson(xr, rank(y)), len(x))))

out.sort(key=lambda t: t[3])
m = len(out)
print(f'\ntested {m} numeric covariates against log2(EBV CN), n={len(x)} donors\n')
print(f'{"covariate":<10} {"pearson":>8} {"spearman":>9} {"p":>10} {"BH q":>10}')
for k, (name, pr, sr, p) in enumerate(out[:15], 1):
    print(f'{name:<10} {pr:>8.3f} {sr:>9.3f} {p:>10.2e} {min(1, p*m/k):>10.2e}')
sig = [t for k, t in enumerate(out, 1) if min(1, t[3] * m / k) < 0.05]
print(f'\n{len(sig)} of {m} covariates significant at BH q<0.05')
if sig:
    print('strongest |spearman|:', max(abs(t[2]) for t in sig).__round__(3))
