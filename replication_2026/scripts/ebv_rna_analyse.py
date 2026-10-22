#!/usr/bin/env python3
"""Comprehensive H1 + H2 + H3 analysis of MAGE EBV RNA-seq quantification."""
import csv, glob, gzip, math, os, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
cn, wt, ir1, ont_cn = {}, {}, {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try: ir1[f[0]] = float(f[5])
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ont_904.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try: ont_cn[f[0]] = float(f[5])
        except ValueError: pass

rows = []
for f in glob.glob(f'{BASE}/rna_mage/*.ebv_rna.tsv'):
    r = list(csv.DictReader(open(f), delimiter='\t'))
    if r and r[0]['donor'] in cn:
        rows.append(r[0])

def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx)**2 for a in x) * sum((b - my)**2 for b in y))
    return num / den if den > 0 else 0.0

def spearman(x, y):
    n = len(x)
    rx = [0]*n; ry = [0]*n
    for idx, (_, i) in enumerate(sorted((v, i) for i, v in enumerate(x))): rx[i] = idx
    for idx, (_, i) in enumerate(sorted((v, i) for i, v in enumerate(y))): ry[i] = idx
    return pearson(rx, ry)

# 1. Multi-replicate donor concordance (technical vs biological thaw)
by_d = defaultdict(list)
for r in rows:
    by_d[r['donor']].append(r)
mult = {d: rs for d, rs in by_d.items() if len(rs) >= 2}
print(f'=== 1. REPLICATE THAW REPRODUCIBILITY ({len(rows)} libraries across {len(by_d)} donors; {len(mult)} multi-replicate donors, {sum(1 for rs in mult.values() if len(rs)==3)} complete 3-run donors) ===')

same_b, diff_b, all_p = [], [], []
for d, rs in mult.items():
    for i in range(len(rs)):
        for j in range(i+1, len(rs)):
            r1, r2 = rs[i], rs[j]
            pair = (float(r1['lytic_late_over_latent_log2']), float(r2['lytic_late_over_latent_log2']),
                    math.log2(float(r1['lytic_ie_rpkm'])+0.1), math.log2(float(r2['lytic_ie_rpkm'])+0.1),
                    math.log2(float(r1['latent_mrna_rpkm'])+0.1), math.log2(float(r2['latent_mrna_rpkm'])+0.1),
                    math.log2(float(r1['lytic_late_rpkm'])+0.1), math.log2(float(r2['lytic_late_rpkm'])+0.1),
                    float(r1['lytic_ie_over_latent_log2']), float(r2['lytic_ie_over_latent_log2']))
            all_p.append(pair)
            if r1['batch'] == r2['batch']:
                same_b.append(pair)
            else:
                diff_b.append(pair)

for label, plist in [('All replicate pairs', all_p), ('Same-batch (technical)', same_b), ('Diff-batch (biological thaw)', diff_b)]:
    if len(plist) >= 3:
        print(f'  {label:28s} (n={len(plist):2d} pairs): r(IE)={pearson([x[2] for x in plist],[x[3] for x in plist]):+.4f}, r(Late)={pearson([x[6] for x in plist],[x[7] for x in plist]):+.4f}, r(IE/Lat)={pearson([x[8] for x in plist],[x[9] for x in plist]):+.4f}, r(Late/Lat)={pearson([x[0] for x in plist],[x[1] for x in plist]):+.4f}, r(Latent)={pearson([x[4] for x in plist],[x[5] for x in plist]):+.4f}')

# 2. Donor-level clean B95-8 H1 correlations (averaging across replicates for multi-run donors)
clean_donors = [d for d in by_d if wt.get(d, 0) <= 0.01]
def d_val(d, col, log_off=None):
    vals = [float(r[col]) for r in by_d[d]]
    m = st.mean(vals)
    return math.log2(m + log_off) if log_off is not None else m

lcn  = [math.log2(cn[d]) for d in clean_donors]
lat  = [d_val(d, 'latent_mrna_rpkm', 0.1) for d in clean_donors]
bart = [d_val(d, 'latent_bart_rpkm', 0.1) for d in clean_donors]
eber = [d_val(d, 'pol3_eber_rpkm', 0.01) for d in clean_donors]
ie   = [d_val(d, 'lytic_ie_rpkm', 0.1) for d in clean_donors]
ear  = [d_val(d, 'lytic_early_rpkm', 0.1) for d in clean_donors]
late = [d_val(d, 'lytic_late_rpkm', 0.1) for d in clean_donors]
rat_late = [d_val(d, 'lytic_late_over_latent_log2') for d in clean_donors]
rat_ie   = [d_val(d, 'lytic_ie_over_latent_log2') for d in clean_donors]
rat_ear  = [d_val(d, 'lytic_early_over_latent_log2') for d in clean_donors]
bz_spl   = [math.log2((st.mean([float(r['bzlf1_spliced_reads']) for r in by_d[d]]) + 1.0) / (st.mean([float(r['ebna_spliced_reads']) for r in by_d[d]]) + 1.0)) for d in clean_donors]
bz_rpm   = [math.log2(st.mean([float(r['bzlf1_spliced_reads'])*1e6/float(r['total_reads']) for r in by_d[d]]) + 0.05) for d in clean_donors]

print(f'\n=== 2. H1 TEST: WGS DNA log2(CN) vs MAGE RNA (n = {len(clean_donors)} unique clean donors) ===')
print(f'  1. log2(CN) vs Latent mRNA RPKM            : r = {pearson(lcn, lat):+.4f} (rho = {spearman(lcn, lat):+.4f})')
print(f'  2. log2(CN) vs BART lncRNA RPKM            : r = {pearson(lcn, bart):+.4f} (rho = {spearman(lcn, bart):+.4f})')
print(f'  3. log2(CN) vs Pol III EBER RPKM           : r = {pearson(lcn, eber):+.4f} (rho = {spearman(lcn, eber):+.4f})')
print(f'  4. log2(CN) vs Lytic IE (BZLF1/BRLF1) RPKM : r = {pearson(lcn, ie):+.4f} (rho = {spearman(lcn, ie):+.4f})')
print(f'  5. log2(CN) vs Spliced BZLF1 (Zta) RPM     : r = {pearson(lcn, bz_rpm):+.4f} (rho = {spearman(lcn, bz_rpm):+.4f})')
print(f'  6. log2(CN) vs Lytic Early RPKM            : r = {pearson(lcn, ear):+.4f} (rho = {spearman(lcn, ear):+.4f})')
print(f'  7. log2(CN) vs Lytic Late Structural RPKM  : r = {pearson(lcn, late):+.4f} (rho = {spearman(lcn, late):+.4f})')
print(f'  8. log2(CN) vs log2(Late / Latent) ratio   : r = {pearson(lcn, rat_late):+.4f} (rho = {spearman(lcn, rat_late):+.4f})')
print(f'  9. log2(CN) vs log2(IE / Latent) ratio     : r = {pearson(lcn, rat_ie):+.4f} (rho = {spearman(lcn, rat_ie):+.4f})')
print(f' 10. log2(CN) vs log2(Early / Latent) ratio  : r = {pearson(lcn, rat_ear):+.4f} (rho = {spearman(lcn, rat_ear):+.4f})')
print(f' 11. log2(CN) vs log2(BZLF1_splice / EBNA_splice): r = {pearson(lcn, bz_spl):+.4f} (rho = {spearman(lcn, bz_spl):+.4f})')

# Within-population demeaned correlation
pop_d = {d: by_d[d][0]['population'] for d in clean_donors}
def demean(vals, groups):
    gm = defaultdict(list)
    for v, g in zip(vals, groups): gm[g].append(v)
    m = {g: st.mean(vs) for g, vs in gm.items()}
    return [v - m[g] for v, g in zip(vals, groups)]

pops = [pop_d[d] for d in clean_donors]
lcn_w = demean(lcn, pops)
print(f'\n=== 3. WITHIN-POPULATION DEMEANED H1 TEST (n = {len(clean_donors)} clean donors across {len(set(pops))} pops) ===')
print(f'  Within-pop log2(CN) vs Latent mRNA RPKM    : r = {pearson(lcn_w, demean(lat, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs Lytic IE RPKM       : r = {pearson(lcn_w, demean(ie, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs Spliced BZLF1 RPM   : r = {pearson(lcn_w, demean(bz_rpm, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs Lytic Early RPKM    : r = {pearson(lcn_w, demean(ear, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs Lytic Late RPKM     : r = {pearson(lcn_w, demean(late, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs log2(IE / Latent)   : r = {pearson(lcn_w, demean(rat_ie, pops)):+.4f}')
print(f'  Within-pop log2(CN) vs log2(Late / Latent) : r = {pearson(lcn_w, demean(rat_late, pops)):+.4f}')

# Also check ONT CN vs MAGE RNA where ONT is available!
ont_d = [d for d in clean_donors if d in ont_cn]
if len(ont_d) >= 15:
    l_ont = [math.log2(ont_cn[d]) for d in ont_d]
    l_30x = [math.log2(cn[d]) for d in ont_d]
    ie_o  = [d_val(d, 'lytic_ie_rpkm', 0.1) for d in ont_d]
    lt_o  = [d_val(d, 'lytic_late_rpkm', 0.1) for d in ont_d]
    la_o  = [d_val(d, 'latent_mrna_rpkm', 0.1) for d in ont_d]
    print(f'\n=== 4. THREE-TIMEPOINT COMPARISON: 2019 30x WGS vs ~2023 ONT vs 2023 MAGE RNA (n = {len(ont_d)} clean donors) ===')
    print(f'  30x log2(CN) vs Lytic IE: r = {pearson(l_30x, ie_o):+.4f} | ONT log2(CN) vs Lytic IE: r = {pearson(l_ont, ie_o):+.4f}')
    print(f'  30x log2(CN) vs Lytic Late: r = {pearson(l_30x, lt_o):+.4f} | ONT log2(CN) vs Lytic Late: r = {pearson(l_ont, lt_o):+.4f}')
    print(f'  30x log2(CN) vs Latent mRNA: r = {pearson(l_30x, la_o):+.4f} | ONT log2(CN) vs Latent mRNA: r = {pearson(l_ont, la_o):+.4f}')

# 5. H3 DIRECT TEST: Does RNA Lytic Reactivation vs Latent mRNA drive the MAGE Host Transcriptome?
Z = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/expression_quants'
with gzip.open(f'{Z}/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz', 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    s_cols = hdr[4:]
    col_idx = [i for i, s in enumerate(s_cols) if s in clean_donors]
    sub_donors = [s_cols[i] for i in col_idx]
    expr_sub = []
    for line in fh:
        p = line.rstrip('\n').split('\t')
        expr_sub.append([float(p[4+i]) for i in col_idx])

x_cn  = [math.log2(cn[d]) for d in sub_donors]
x_ie  = [d_val(d, 'lytic_ie_rpkm', 0.1) for d in sub_donors]
x_lt  = [d_val(d, 'lytic_late_rpkm', 0.1) for d in sub_donors]
x_rat = [d_val(d, 'lytic_ie_over_latent_log2') for d in sub_donors]
x_lat = [d_val(d, 'latent_mrna_rpkm', 0.1) for d in sub_donors]

r_cn_g  = [pearson(x_cn, g) for g in expr_sub]
r_ie_g  = [pearson(x_ie, g) for g in expr_sub]
r_lt_g  = [pearson(x_lt, g) for g in expr_sub]
r_rat_g = [pearson(x_rat, g) for g in expr_sub]
r_lat_g = [pearson(x_lat, g) for g in expr_sub]

print(f'\n=== 5. H3 HOST TRANSCRIPTOME VECTOR ALIGNMENT across 20,154 Host Genes (n = {len(sub_donors)} clean donors) ===')
print(f'  Corr across 20,154 genes of r(Gene, WGS_CN) with r(Gene, RNA_Lytic_IE)     : r = {pearson(r_cn_g, r_ie_g):+.4f}')
print(f'  Corr across 20,154 genes of r(Gene, WGS_CN) with r(Gene, RNA_Lytic_Late)   : r = {pearson(r_cn_g, r_lt_g):+.4f}')
print(f'  Corr across 20,154 genes of r(Gene, WGS_CN) with r(Gene, log2_IE_over_Lat) : r = {pearson(r_cn_g, r_rat_g):+.4f}')
print(f'  Corr across 20,154 genes of r(Gene, WGS_CN) with r(Gene, RNA_Latent_mRNA)  : r = {pearson(r_cn_g, r_lat_g):+.4f}')
print(f'  Mean |r| across 20,154 host genes: WGS_CN={st.mean([abs(v) for v in r_cn_g]):.4f}, Lytic_IE={st.mean([abs(v) for v in r_ie_g]):.4f}, Lytic_Late={st.mean([abs(v) for v in r_lt_g]):.4f}, Latent_mRNA={st.mean([abs(v) for v in r_lat_g]):.4f}')
