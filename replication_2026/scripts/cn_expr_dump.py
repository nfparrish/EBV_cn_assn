"""Dump per-gene CN-expression correlations (marginal, within-population, and
permuted nulls for each) so the figure can show the full distribution, not just a count.
Reproduces cn_expression.py and within_pop.py in one pass.
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

def prep(v):
    m = st.mean(v); d = [p-m for p in v]
    return d, math.sqrt(sum(p*p for p in d))

out = open('/work/nfp8/EBV/agg/cn_expression_per_gene.tsv','w')
out.write('gene\trho_marginal\trho_marginal_perm\trho_withinpop\trho_withinpop_perm\n')

topdump = open('/work/nfp8/EBV/agg/cn_expression_topgenes.tsv','w')

with gzip.open(BED,'rt') as fh:
    hdr = next(fh).rstrip('\n').split('\t')
    sid = hdr[4:]
    cols = [j for j,s in enumerate(sid) if s in cn]
    x = [cn[sid[j]] for j in cols]
    n = len(x)
    groups = defaultdict(list)
    for j,s in enumerate(sid):
        if s in cn and s in pop: groups[pop[s]].append(j)
    groups = {k:v for k,v in groups.items() if len(v) >= 20}
    n_eff = sum(len(v) for v in groups.values())
    print(f'marginal donors: {n}; within-pop populations: {len(groups)}; within-pop donors: {n_eff}', flush=True)

    xr = rank(x); dx, sx = prep(xr)
    random.seed(0)
    xp = xr[:]; random.shuffle(xp); dxp, sxp = prep(xp)

    gx, gxp = {}, {}
    for g, js in groups.items():
        r_ = rank([cn[sid[j]] for j in js]); gx[g] = prep(r_)
        pm = r_[:]; random.shuffle(pm); gxp[g] = prep(pm)

    thr_m = 3.5/math.sqrt(n); thr_w = 3.5/math.sqrt(n_eff)
    print(f'thr_marginal={thr_m:.4f} thr_withinpop={thr_w:.4f}', flush=True)

    rows = []
    ng = 0
    for line in fh:
        p = line.rstrip('\n').split('\t'); e = p[4:]
        try: y = [float(e[j]) for j in cols]
        except (ValueError, IndexError): continue
        dy, sy = prep(rank(y))
        if not sy: continue
        rm  = sum(a*b for a,b in zip(dx,dy))/(sx*sy)
        rmp = sum(a*b for a,b in zip(dxp,dy))/(sxp*sy)
        zs = zps = w = 0.0; ok = True
        for g, js in groups.items():
            try: yg = [float(e[j]) for j in js]
            except (ValueError, IndexError): ok=False; break
            dyg, syg = prep(rank(yg))
            if not syg: ok=False; break
            d1,s1 = gx[g]; d2,s2 = gxp[g]
            r  = sum(a*b for a,b in zip(d1,dyg))/(s1*syg)
            rp = sum(a*b for a,b in zip(d2,dyg))/(s2*syg)
            k = len(js)-3
            r = max(min(r,.999),-.999); rp = max(min(rp,.999),-.999)
            zs += k*0.5*math.log((1+r)/(1-r)); zps += k*0.5*math.log((1+rp)/(1-rp)); w += k
        if not ok or w<=0: continue
        rw = math.tanh(zs/w); rwp = math.tanh(zps/w)
        ng += 1
        out.write(f'{p[3]}\t{rm:.5f}\t{rmp:.5f}\t{rw:.5f}\t{rwp:.5f}\n')
        rows.append((abs(rm), p[3], y))
        if ng % 2000 == 0: print(f'  {ng} genes', flush=True)
out.close()
rows.sort(reverse=True)
# dump expression vectors for the top 3 genes plus the donor CN, for a scatter panel
topdump.write('sample\tlog2_cn\t' + '\t'.join(g for _,g,_ in rows[:3]) + '\n')
for i,j in enumerate(cols):
    topdump.write(f'{sid[j]}\t{x[i]:.5f}\t' + '\t'.join(f'{rows[k][2][i]:.5f}' for k in range(3)) + '\n')
topdump.close()
print(f'genes written: {ng}')
