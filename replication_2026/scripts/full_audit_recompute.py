import csv, glob, math, os, random, statistics as st
from collections import defaultdict, Counter

BASE = '/work/nfp8/EBV'

pop, spop, sex_ped, fam_id = {}, {}, {}, {}
trios = []
for line in open(f'{BASE}/metadata/ped3202.txt'):
    f = line.split()
    if len(f) >= 7 and f[1] != 'SampleID':
        s = f[1]
        fam_id[s] = f[0]
        sex_ped[s] = f[4]
        pop[s] = f[5]    # 26 populations
        spop[s] = f[6]   # 5 superpopulations
        if f[2] != '0' and f[3] != '0':
            trios.append((s, f[2], f[3]))

mage_meta = {}
for m in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'), delimiter='\t'):
    mage_meta[m['sample_kgpID']] = m

cn, wt, auto_d, ebv_d = {}, {}, {}, {}
with open(f'{BASE}/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            v = float(r['ebv_cn_depth_ratio']); w = float(r['wt_fraction'])
            if v > 0:
                s = r['sample']
                cn[s] = v; wt[s] = w
                auto_d[s] = float(r['auto_depth']); ebv_d[s] = float(r['ebv_depth'])
        except ValueError: pass

unrel_2502 = set(cn.keys())
for f in os.listdir(f'{BASE}/cn'):
    if not f.endswith('.tsv'): continue
    s = f[:-4]
    if s in cn: continue
    rows = list(csv.DictReader(open(f'{BASE}/cn/{f}'), delimiter='\t'))
    if rows:
        try:
            v = float(rows[0]['ebv_cn_depth_ratio']); w = float(rows[0]['wt_fraction'])
            if v > 0:
                cn[s] = v; wt[s] = w
                auto_d[s] = float(rows[0]['auto_depth']); ebv_d[s] = float(rows[0]['ebv_depth'])
        except ValueError: pass

coinf = {s for s, v in wt.items() if v > 0.01}
clean_2502 = {s for s in unrel_2502 if s not in coinf}
clean_all  = {s for s in cn if s not in coinf}

print(f'=== 1. COHORT COMPLETENESS ===')
print(f'30x unrelated CN: {len(unrel_2502)}/2504 (clean: {len(clean_2502)}, co-infected wt>1%: {len(unrel_2502 & coinf)})')
print(f'30x total CN    : {len(cn)}/3202 (clean: {len(clean_all)}, co-infected wt>1%: {len(coinf)})')

RU = 22737 / 3072.0
ir1, uniq_d = {}, {}
for f in glob.glob(f'{BASE}/ir1/*.tsv'):
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    if rows:
        try:
            s = rows[0]['sample']
            v = float(rows[0]['ir1_copies_vs_unique'])
            u = float(rows[0]['unique_depth'])
            if v > 0: ir1[s] = v * RU
            if u > 0: uniq_d[s] = u
        except (ValueError, KeyError): pass

cn_uniq = {s: 2.0 * uniq_d[s] / auto_d[s] for s in uniq_d if s in auto_d and auto_d[s] > 0}
ratios_u = [cn_uniq[s] / cn[s] for s in clean_2502 if s in cn_uniq]
print(f'\n=== 2. UNIQUE-SEQUENCE CN vs STANDARD CN (n={len(ratios_u)} clean lines) ===')
print(f'  Standard CN (clean 2,449): min {min(cn[s] for s in clean_2502):.2f}, median {st.median(cn[s] for s in clean_2502):.2f}, max {max(cn[s] for s in clean_2502):.2f}')
print(f'  Unique-seq CN (clean)    : min {min(cn_uniq[s] for s in clean_2502):.2f}, median {st.median(cn_uniq[s] for s in clean_2502):.2f}, max {max(cn_uniq[s] for s in clean_2502):.2f}')
print(f'  Ratio cn_uniq / cn_std   : median {st.median(ratios_u):.4f} (inv={1/st.median(ratios_u):.4f}), mean {st.mean(ratios_u):.4f}, sd {st.pstdev(ratios_u):.4f}')

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    da = [x - ma for x in a]; db = [y - mb for y in b]
    den = math.sqrt(sum(x*x for x in da) * sum(y*y for y in db))
    return sum(x*y for x, y in zip(da, db)) / den if den else 0.0

def rank(vals):
    o = sorted(range(len(vals)), key=lambda i: vals[i]); r = [0.0]*len(vals); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and vals[o[j+1]] == vals[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r

def spearman(a, b):
    return pearson(rank(a), rank(b))

def kw_eta(groups):
    flat = [v for g in groups.values() for v in g]; n = len(flat)
    rk = rank(flat); pos = 0; H = 0.0
    for g in groups.values():
        rs = sum(rk[pos:pos+len(g)]); pos += len(g)
        H += rs * rs / len(g)
    H = 12.0 / (n * (n + 1)) * H - 3 * (n + 1); k = len(groups)
    return max(0.0, (H - k + 1) / (n - k)), H, k, n

def anova_r2(groups):
    flat = [v for g in groups.values() for v in g]; gm = st.mean(flat)
    sst = sum((v - gm)**2 for v in flat)
    ssb = sum(len(g) * (st.mean(g) - gm)**2 for g in groups.values())
    return ssb / sst if sst else 0.0

s_u = sorted(clean_2502 & set(cn_uniq))
print(f'  Pearson r (log2): {pearson([math.log2(cn[s]) for s in s_u], [math.log2(cn_uniq[s]) for s in s_u]):.5f}')

print(f'\n=== 3. CN POPULATION STRUCTURE: ALL 2,449 CLEAN 1kGP vs 715 CLEAN MAGE ===')
for label, subset in [('All 2,449 clean 1kGP (ped3202)', clean_2502), ('MAGE 715 clean subset', clean_2502 & set(mage_meta))]:
    g_pop = defaultdict(list); g_cont = defaultdict(list)
    for s in subset:
        if s in pop:
            g_pop[pop[s]].append(math.log2(cn[s]))
            g_cont[spop[s]].append(math.log2(cn[s]))
    eta_p, H_p, k_p, n_p = kw_eta(g_pop)
    eta_c, H_c, k_c, n_c = kw_eta(g_cont)
    r2_p = anova_r2(g_pop); r2_c = anova_r2(g_cont)
    print(f'  {label} (n={n_p}, {k_p} pops):')
    print(f'    Population  (26): KW eta^2 = {eta_p:.4f}, ANOVA R^2 = {r2_p:.4f}')
    print(f'    Continental  (5): KW eta^2 = {eta_c:.4f}, ANOVA R^2 = {r2_c:.4f}')
    meds = sorted(((2**st.median(v), p, len(v)) for p, v in g_pop.items()), reverse=True)
    print(f'    Top 4 pops   : ' + ', '.join(f'{p} {m:.1f} (n={n})' for m, p, n in meds[:4]))
    print(f'    Bottom 4 pops: ' + ', '.join(f'{p} {m:.1f} (n={n})' for m, p, n in meds[-4:]))
    eur_meds = sorted(((2**st.median(g_pop[p]), p, len(g_pop[p])) for p in g_pop if spop[next(s for s in subset if pop[s]==p)]=='EUR'), reverse=True)
    print(f'    EUR pops     : ' + ', '.join(f'{p} {m:.1f} (n={n})' for m, p, n in eur_meds))
