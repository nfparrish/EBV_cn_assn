"""Is the CN-expression association independent of ancestry?

CN is population-structured (eta-sq 0.126), so a genome-wide CN-expression correlation
could reflect ancestry driving both. Control: rank-correlate expression against CN
*within* each of the 26 populations, then average the within-population correlations
(Fisher z). Ancestry is approximately constant within a population, so a surviving
association is CN acting independently of ancestry.
Permuted null: shuffle CN within each population, preserving group structure.
"""
import csv, glob, gzip, math, random, statistics as st
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

pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'),
    delimiter='\t')}

BED = ('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/'
       'expression_quants/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz')

def rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r=[0.0]*len(v); i=0
    while i < len(o):
        j=i
        while j+1<len(o) and v[o[j+1]]==v[o[i]]: j+=1
        for k in range(i,j+1): r[o[k]]=(i+j)/2+1
        i=j+1
    return r

with gzip.open(BED, 'rt') as fh:
    sid = next(fh).rstrip('\n').split('\t')[4:]
    groups = defaultdict(list)
    for j, s in enumerate(sid):
        if s in cn and s in pop: groups[pop[s]].append(j)
    groups = {k: v for k, v in groups.items() if len(v) >= 20}
    print(f'{len(groups)} populations, {sum(len(v) for v in groups.values())} donors '
          f'(median {int(st.median([len(v) for v in groups.values()]))} per population)')

    random.seed(0)
    gx, gxp = {}, {}
    for g, js in groups.items():
        xr = rank([cn[sid[j]] for j in js]); gx[g] = xr
        pm = xr[:]; random.shuffle(pm); gxp[g] = pm

    def prep(v):
        m = st.mean(v); d = [p-m for p in v]
        return d, math.sqrt(sum(p*p for p in d))

    pg  = {g: prep(v) for g, v in gx.items()}
    pgp = {g: prep(v) for g, v in gxp.items()}

    n_eff = sum(len(v) for v in groups.values())
    thr = 3.5/math.sqrt(n_eff)
    hits = ph = ngene = 0; top = []
    for line in fh:
        p = line.rstrip('\n').split('\t'); e = p[4:]
        zs = zps = w = 0.0; ok = True
        for g, js in groups.items():
            try: y = [float(e[j]) for j in js]
            except (ValueError, IndexError): ok = False; break
            dy, sy = prep(rank(y))
            if not sy: ok = False; break
            dx, sx = pg[g]; dxp, sxp = pgp[g]
            r  = sum(a*b for a, b in zip(dx, dy))/(sx*sy)
            rp = sum(a*b for a, b in zip(dxp, dy))/(sxp*sy)
            k = len(js)-3
            r = max(min(r, .999), -.999); rp = max(min(rp, .999), -.999)
            zs += k*0.5*math.log((1+r)/(1-r)); zps += k*0.5*math.log((1+rp)/(1-rp)); w += k
        if not ok or w <= 0: continue
        ngene += 1
        rm  = math.tanh(zs/w); rpm = math.tanh(zps/w)
        if abs(rm)  > thr: hits += 1
        if abs(rpm) > thr: ph += 1
        top.append((abs(rm), rm, p[3]))

top.sort(reverse=True)
print(f'\ngenes tested: {ngene}   |rho| threshold: {thr:.3f}')
print(f'  within-population hits : {hits}  ({100*hits/ngene:.1f}%)')
print(f'  permuted-within-pop    : {ph}  ({100*ph/ngene:.1f}%)   <- null')
print(f'  enrichment             : {hits/max(ph,1):.1f}x\n')
print('strongest within-population (ancestry-controlled) associations:')
for a, r, g in top[:10]:
    print(f'  {g:<22} mean within-pop rho = {r:+.3f}')
