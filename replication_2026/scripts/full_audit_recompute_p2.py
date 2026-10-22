import csv, glob, gzip, math, os, random, statistics as st
from collections import defaultdict, Counter

BASE = '/work/nfp8/EBV'
RU = 22737 / 3072.0

pop, spop, fam_id = {}, {}, {}
trios = []
for line in open(f'{BASE}/metadata/ped3202.txt'):
    f = line.split()
    if len(f) >= 7 and f[1] != 'SampleID':
        s = f[1]; fam_id[s] = f[0]; pop[s] = f[5]; spop[s] = f[6]
        if f[2] != '0' and f[3] != '0':
            trios.append((s, f[2], f[3]))

cn, wt = {}, {}
for f in os.listdir(f'{BASE}/cn'):
    if not f.endswith('.tsv'): continue
    rows = list(csv.DictReader(open(f'{BASE}/cn/{f}'), delimiter='\t'))
    if rows:
        try:
            v = float(rows[0]['ebv_cn_depth_ratio']); w = float(rows[0]['wt_fraction'])
            if v > 0: cn[rows[0]['sample']] = v; wt[rows[0]['sample']] = w
        except ValueError: pass

unrel_2502 = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t') if r['sample'] in cn}
coinf = {s for s, v in wt.items() if v > 0.01}
clean_2502 = unrel_2502 - coinf

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

def kw_eta(groups):
    flat = [v for g in groups.values() for v in g]; n = len(flat)
    rk = rank(flat); pos = 0; H = 0.0
    for g in groups.values():
        rs = sum(rk[pos:pos+len(g)]); pos += len(g)
        H += rs * rs / len(g)
    H = 12.0 / (n * (n + 1)) * H - 3 * (n + 1); k = len(groups)
    return max(0.0, (H - k + 1) / (n - k)), H, k, n

def cohen(a, b):
    p = math.sqrt((st.pvariance(a)*len(a) + st.pvariance(b)*len(b)) / (len(a) + len(b)))
    return (st.mean(a) - st.mean(b)) / p if p else 0.0

# === 4. HERITABILITY FULL BREAKDOWN ===
full_trios = [t for t in trios if all(m in cn for m in t) and not any(m in coinf for m in t)]
print(f'=== 4. HERITABILITY & SPOUSE CONTROL (n={len(full_trios)} clean trios) ===')
unique_spouses = list({(f, m) for c, f, m in full_trios})
print(f'  Unique father-mother spouse pairs: {len(unique_spouses)} (out of {len(full_trios)} trios -> {len(full_trios)-len(unique_spouses)} sibling repeats)')
for tname in ('raw', 'log2', 'rank'):
    if tname == 'raw': vals = dict(cn)
    elif tname == 'log2': vals = {s: math.log2(v) for s, v in cn.items()}
    else:
        ks = sorted(cn); rks = rank([cn[k] for k in ks]); vals = {k: r/len(ks) for k, r in zip(ks, rks)}
    # Pooled
    r_mid = pearson([vals[c] for c, f, m in full_trios], [(vals[f]+vals[m])/2 for c, f, m in full_trios])
    r_sp  = pearson([vals[c] for c, f, m in full_trios]*2, [vals[f] for c, f, m in full_trios]+[vals[m] for c, f, m in full_trios])
    r_pp  = pearson([vals[f] for c, f, m in full_trios], [vals[m] for c, f, m in full_trios])
    r_pp_uniq = pearson([vals[f] for f, m in unique_spouses], [vals[m] for f, m in unique_spouses])
    # Within-population (demeaned by actual population f[5])
    trio_pops = set(pop[c] for c, f, m in full_trios)
    pmeans = {p: st.mean([vals[s] for s in cn if pop.get(s)==p and s not in coinf]) for p in trio_pops}
    wvals = {s: vals[s] - pmeans[pop[s]] for t in full_trios for s in t}
    r_mid_w = pearson([wvals[c] for c, f, m in full_trios], [(wvals[f]+wvals[m])/2 for c, f, m in full_trios])
    r_sp_w  = pearson([wvals[c] for c, f, m in full_trios]*2, [wvals[f] for c, f, m in full_trios]+[wvals[m] for c, f, m in full_trios])
    r_pp_w  = pearson([wvals[f] for c, f, m in full_trios], [wvals[m] for c, f, m in full_trios])
    r_pp_w_uniq = pearson([wvals[f] for f, m in unique_spouses], [wvals[m] for f, m in unique_spouses])
    print(f'  Transform={tname:<4}:')
    print(f'    Pooled     : r_mid={r_mid:+.4f}, r_sp={r_sp:+.4f}, r_pp={r_pp:+.4f} (uniq spouses: {r_pp_uniq:+.4f}), 2(PO-PP)={2*(r_sp-r_pp):+.4f}')
    print(f'    Within-pop : r_mid={r_mid_w:+.4f}, r_sp={r_sp_w:+.4f}, r_pp={r_pp_w:+.4f} (uniq spouses: {r_pp_w_uniq:+.4f}), 2(PO-PP)={2*(r_sp_w-r_pp_w):+.4f}')

# === 5. IR1 REPEAT COPY NUMBER ON ALL 2,449 CLEAN 1kGP + GC CONTROL (n=2,449) + ONT IR1 (n=907) ===
ir1 = {}
for f in glob.glob(f'{BASE}/ir1/*.tsv'):
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    if rows:
        try:
            v = float(rows[0]['ir1_copies_vs_unique'])
            if v > 0: ir1[rows[0]['sample']] = v * RU
        except (ValueError, KeyError): pass

gcb = {}
for f in glob.glob(f'{BASE}/gcbias/*.tsv'):
    p = open(f).read().rstrip('\n').split('\t')
    if len(p) >= 5 and p[4] != 'NA':
        try: gcb[p[0]] = float(p[4])
        except ValueError: pass

ont_ir1, ont_span_med, ont_span_max = {}, {}, {}
for f in glob.glob(f'{BASE}/ont_ir1/*.tsv'):
    rows = list(csv.DictReader(open(f), delimiter='\t'))
    if rows:
        s = rows[0]['sample']
        try:
            v = float(rows[0]['units_by_depth'])
            if v > 0: ont_ir1[s] = v
        except (ValueError, KeyError): pass
        try:
            ont_span_med[s] = float(rows[0]['median_units_spanned'])
            ont_span_max[s] = float(rows[0]['max_units_spanned'])
        except (ValueError, KeyError): pass

clean_ir1 = {s: v for s, v in ir1.items() if s in clean_2502 and s in pop}
print(f'\n=== 5. IR1 REPEAT COPY NUMBER (n={len(clean_ir1)} clean 1kGP, ONT n={len(ont_ir1)}) ===')
g_pop_ir1 = defaultdict(list); g_cont_ir1 = defaultdict(list)
for s, v in clean_ir1.items():
    g_pop_ir1[pop[s]].append(v)
    g_cont_ir1[spop[s]].append(v)
eta_p_ir1, _, _, _ = kw_eta(g_pop_ir1)
eta_c_ir1, _, _, _ = kw_eta(g_cont_ir1)
print(f'  Overall clean IR1: median {st.median(clean_ir1.values()):.2f}u, mean {st.mean(clean_ir1.values()):.2f}u, sd {st.pstdev(clean_ir1.values()):.2f}u, range {min(clean_ir1.values()):.2f}-{max(clean_ir1.values()):.2f}')
print(f'  KW eta^2: Population (26) = {eta_p_ir1:.4f}, Continental (5) = {eta_c_ir1:.4f} (ratio cont/pop = {eta_c_ir1/eta_p_ir1:.3f})')
for p in ['CEU', 'FIN', 'IBS', 'GBR', 'TSI']:
    print(f'    EUR panel {p:<4}: n={len(g_pop_ir1[p]):<4} median={st.median(g_pop_ir1[p]):.2f}u, sd={st.pstdev(g_pop_ir1[p]):.2f}u')
for c in ['EUR', 'EAS', 'AMR', 'SAS', 'AFR']:
    print(f'    Continent {c:<4}: n={len(g_cont_ir1[c]):<4} median={st.median(g_cont_ir1[c]):.2f}u')

# GC residualisation on ALL 2,449 clean lines (not 182!)
ss_gc = sorted(s for s in clean_ir1 if s in gcb)
y_ir = [clean_ir1[s] for s in ss_gc]; x_gc = [gcb[s] for s in ss_gc]
mx, my = st.mean(x_gc), st.mean(y_ir)
b_gc = sum((a-mx)*(c-my) for a, c in zip(x_gc, y_ir)) / sum((a-mx)**2 for a in x_gc)
res_gc = [c - (my + b_gc*(a-mx)) for a, c in zip(x_gc, y_ir)]
ve_gc = (b_gc**2) * st.pvariance(x_gc) / st.pvariance(y_ir)
eur4 = {'CEU', 'FIN', 'IBS', 'GBR'}
raw_e4 = [clean_ir1[s] for s in ss_gc if pop[s] in eur4]; raw_o4 = [clean_ir1[s] for s in ss_gc if pop[s] not in eur4]
adj_e4 = [r for s, r in zip(ss_gc, res_gc) if pop[s] in eur4]; adj_o4 = [r for s, r in zip(ss_gc, res_gc) if pop[s] not in eur4]
raw_eu = [clean_ir1[s] for s in ss_gc if spop[s]=='EUR']; raw_oe = [clean_ir1[s] for s in ss_gc if spop[s]!='EUR']
adj_eu = [r for s, r in zip(ss_gc, res_gc) if spop[s]=='EUR']; adj_oe = [r for s, r in zip(ss_gc, res_gc) if spop[s]!='EUR']
print(f'  GC control on ALL {len(ss_gc)} clean lines: GC explains {100*ve_gc:.2f}% of IR1 variance')
print(f'    EUR (all 5 incl TSI) vs non-EUR : raw d={cohen(raw_eu, raw_oe):+.3f} -> GC-adj d={cohen(adj_eu, adj_oe):+.3f} ({100*abs(cohen(adj_eu, adj_oe))/abs(cohen(raw_eu, raw_oe)):.1f}% survives)')
print(f'    4 low-EUR (CEU/FIN/GBR/IBS) vs rest: raw d={cohen(raw_e4, raw_o4):+.3f} -> GC-adj d={cohen(adj_e4, adj_o4):+.3f} ({100*abs(cohen(adj_e4, adj_o4))/abs(cohen(raw_e4, raw_o4)):.1f}% survives)')

# ONT IR1 on ALL 907 lines (not 40 or 705!)
ss_ont = sorted(s for s in ont_ir1 if s in clean_ir1)
r_ont_ir1 = pearson([clean_ir1[s] for s in ss_ont], [ont_ir1[s] for s in ss_ont])
ont_eu = [ont_ir1[s] for s in ss_ont if spop[s]=='EUR']; ont_oe = [ont_ir1[s] for s in ss_ont if spop[s]!='EUR']
ont_e4 = [ont_ir1[s] for s in ss_ont if pop[s] in eur4]; ont_o4 = [ont_ir1[s] for s in ss_ont if pop[s] not in eur4]
print(f'  ONT IR1 (n={len(ss_ont)} clean paired with 30x): r = {r_ont_ir1:.4f}, median ratio ONT/30x = {st.median([ont_ir1[s]/clean_ir1[s] for s in ss_ont]):.4f}')
print(f'    ONT EUR vs non-EUR: d = {cohen(ont_eu, ont_oe):+.3f}; 4 low-EUR vs rest: d = {cohen(ont_e4, ont_o4):+.3f}')
for p in ['CEU', 'FIN', 'IBS', 'GBR', 'TSI']:
    vals_p = [ont_ir1[s] for s in ss_ont if pop[s]==p]
    if vals_p: print(f'    ONT {p:<4}: n={len(vals_p):<3} median={st.median(vals_p):.2f}u')
print(f'    ONT structural span cap check: max(median_units_spanned)={max(ont_span_med.values()):.2f}, max(max_units_spanned)={max(ont_span_max.values()):.2f} (reference=7.40)')
