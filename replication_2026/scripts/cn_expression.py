"""Does EBV copy number move host gene expression in MAGE?

Closes the confounding argument. MAGE's expression ANOVA regresses out only sex+batch,
then attributes residual variance to population labels. CN structures with population
(eta-sq 0.126) and is NOT absorbed by batch (0.003, n.s.). If CN also associates with
host expression, then variance MAGE attributes to ancestry is partly viral.

Null calibration by permuting the CN vector, so the hit count is interpretable.
"""
import csv, glob, gzip, math, random, statistics as st

cn = {}
for f in glob.glob('/work/nfp8/EBV/cn/*.tsv'):
    with open(f) as fh:
        r = list(csv.DictReader(fh, delimiter='\t'))
    if r:
        try:
            v = float(r[0]['ebv_cn_depth_ratio'])
            if v > 0:
                cn[r[0]['sample']] = math.log2(v)
        except ValueError:
            pass

BED = ('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/'
       'expression_quants/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz')
with gzip.open(BED, 'rt') as fh:
    hdr = next(fh).rstrip('\n').split('\t')
    sid = hdr[4:]
    keep = [(j, s) for j, s in enumerate(sid) if s in cn]
    cols = [j for j, _ in keep]
    x = [cn[s] for _, s in keep]
    print(f'{len(x)} donors with both CN and expression; {len(sid)} in matrix')

    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(o):
            j = i
            while j+1 < len(o) and v[o[j+1]] == v[o[i]]: j += 1
            for k in range(i, j+1): r[o[k]] = (i+j)/2+1
            i = j+1
        return r

    xr = rank(x); n = len(x)
    mx = st.mean(xr); dx = [p-mx for p in xr]; sx = math.sqrt(sum(p*p for p in dx))
    random.seed(0)
    xperm = xr[:]; random.shuffle(xperm)
    mp = st.mean(xperm); dp = [p-mp for p in xperm]; sp_ = math.sqrt(sum(p*p for p in dp))

    def corr(dv, sv, d2, s2):
        return sum(a*b for a, b in zip(dv, d2))/(sv*s2) if sv and s2 else 0.0

    thr = 3.5 / math.sqrt(n)        # |rho| giving roughly p<5e-4
    hits = perm_hits = ngene = 0
    top = []
    for line in fh:
        p = line.rstrip('\n').split('\t')
        try:
            y = [float(p[4:][j]) for j in cols]
        except (ValueError, IndexError):
            continue
        yr = rank(y); my = st.mean(yr); dy = [q-my for q in yr]
        sy = math.sqrt(sum(q*q for q in dy))
        r = corr(dx, sx, dy, sy); rp = corr(dp, sp_, dy, sy)
        ngene += 1
        if abs(r) > thr: hits += 1
        if abs(rp) > thr: perm_hits += 1
        top.append((abs(r), r, p[3]))

top.sort(reverse=True)
print(f'\ngenes tested: {ngene}   |rho| threshold: {thr:.3f} (~p<5e-4)')
print(f'  observed hits : {hits}  ({100*hits/ngene:.1f}% of genes)')
print(f'  permuted hits : {perm_hits}  ({100*perm_hits/ngene:.1f}%)  <- null calibration')
print(f'  enrichment    : {hits/max(perm_hits,1):.1f}x\n')
print('strongest associations with log2 EBV CN:')
for a, r, g in top[:12]:
    print(f'  {g:<22} rho = {r:+.3f}')
