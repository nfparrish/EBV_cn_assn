#!/usr/bin/env python3
"""Aim 2 (Designs B, C, D): Viral genetic influence on EBV copy number and viral gene expression:
  - Design B: Within-population Viral CN-QTLs (717 MAGE + 1,732 non-MAGE replication = 2,449 clean 1kGP)
              & Within-population+batch trans-lytic and cis-gene Viral eQTLs in MAGE (n=717).
  - Design C: IR1 (BamHI-W) repeat copy number as a structural cis-eQTL (W1-W2, EBNA splices, BHRF1, EBNA-2).
  - Design D: Rare viral regulatory variant burden (Zp/Rp, oriP, oriLyt) & Wild-type LF2 co-infection effect on Late/IE progression.
"""
import csv, math, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
ped = {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    if line.startswith('FamilyID'): continue
    f = line.split()
    if len(f) >= 7: ped[f[1]] = (f[5], f[6]) # pop, superpop

cn, wt, ir1, ont_ir1 = {}, {}, {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5:
        try: ir1[f[0]] = float(f[4]) * (22737.0/3072.0)
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ont_ir1_907.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5:
        try: ont_ir1[f[0]] = float(f[3])
        except ValueError: pass

unrel = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t') if r['sample'] in cn}
clean_wgs = sorted([s for s in unrel if wt[s] <= 0.01 and s in ped])

rna_rows = [r for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_rna_779.tsv'), delimiter='\t') if r['donor'] in cn]
by_d = defaultdict(list)
for r in rna_rows: by_d[r['donor']].append(r)
clean_mage = sorted([d for d in by_d if wt.get(d, 0) <= 0.01 and d in ped])
coinf_mage = sorted([d for d in by_d if wt.get(d, 0) > 0.01 and d in ped])
non_mage_wgs = sorted([s for s in clean_wgs if s not in by_d])

# Load per-gene RPKM for MAGE donors
run_gene = defaultdict(dict)
gene_class = {}
for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_genes_779.tsv'), delimiter='\t'):
    run_gene[r['run']][r['gene']] = float(r['rpkm'])
    gene_class[r['gene']] = r['kinetic_class']

donor_gene = {d: {g: st.mean([run_gene[r['run']][g] for r in by_d[d]]) for g in gene_class} for d in by_d}

def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx)**2 for a in x) * sum((b - my)**2 for b in y))
    return num / den if den > 0 else 0.0

def pv(r, n):
    if n <= 3 or abs(r) >= 1.0: return 0.0
    t = abs(r) * math.sqrt((n - 2) / (1.0 - r*r))
    return math.erfc(t / math.sqrt(2.0))

def demean(vals, groups):
    gm = defaultdict(list)
    for v, g in zip(vals, groups): gm[g].append(v)
    m = {g: st.mean(vs) for g, vs in gm.items()}
    return [v - m[g] for v, g in zip(vals, groups)]

# ==============================================================================
# DESIGN C: IR1 (BamHI-W) Repeat Copy Number as a Structural cis-QTL
# ==============================================================================
m_ir1 = [d for d in clean_mage if d in ir1]
pops_m = [ped[d][0] for d in m_ir1]
bats_m = [by_d[d][0]['batch'] for d in m_ir1]
x_ir1_raw = [ir1[d] for d in m_ir1]
x_ir1_wb  = demean(demean(x_ir1_raw, pops_m), bats_m)

def d_metric(d, fn):
    return st.mean([fn(r) for r in by_d[d]])

ir1_rpm     = [math.log2(d_metric(d, lambda r: float(r['ir1_mapq0_reads'])*1e6/float(r['total_reads'])) + 0.1) for d in m_ir1]
ir1_ov_lat  = [math.log2((d_metric(d, lambda r: float(r['ir1_mapq0_reads'])*1e6/float(r['total_reads'])) + 0.1) / (d_metric(d, lambda r: float(r['latent_mrna_rpkm'])) + 0.1)) for d in m_ir1]
ebna_spl    = [math.log2(d_metric(d, lambda r: float(r['ebna_spliced_reads'])*1e6/float(r['total_reads'])) + 0.1) for d in m_ir1]
ebna_ov_lat = [math.log2((d_metric(d, lambda r: float(r['ebna_spliced_reads'])) + 1.0) / (d_metric(d, lambda r: float(r['latent_mrna_reads'])) + 1.0)) for d in m_ir1]
ebna2_rpkm  = [math.log2(donor_gene[d]['EBNA-2'] + 0.1) for d in m_ir1]
bhrf1_rpkm  = [math.log2(donor_gene[d]['BHRF1'] + 0.1) for d in m_ir1]
bhrf1_ov_ie = [math.log2((donor_gene[d]['BHRF1'] + 0.1) / (d_metric(d, lambda r: float(r['lytic_ie_rpkm'])) + 0.1)) for d in m_ir1]
lmp1_rpkm   = [math.log2(donor_gene[d]['LMP-1'] + 0.1) for d in m_ir1]
ie_lat_rat  = [d_metric(d, lambda r: float(r['lytic_ie_over_latent_log2'])) for d in m_ir1]
lcn_m       = [math.log2(cn[d]) for d in m_ir1]

print(f'=== AIM 2.3 (DESIGN C): IR1 (BamHI-W) REPEAT COPY NUMBER AS A STRUCTURAL CIS-QTL (n = {len(m_ir1)} clean MAGE donors) ===')
for name, vec in [
    ('IR1 (W1-W2) Repeat RNA RPM', ir1_rpm),
    ('log2(IR1 W1-W2 RNA / Latent mRNA)', ir1_ov_lat),
    ('EBNA W1-W2/Cp/Wp Spliced Reads RPM', ebna_spl),
    ('log2(EBNA_spliced / Latent_mRNA_reads)', ebna_ov_lat),
    ('BHRF1 (W-spliced + lytic vBCL2) RPKM', bhrf1_rpkm),
    ('log2(BHRF1 / Lytic IE) [Latent W-BHRF1]', bhrf1_ov_ie),
    ('EBNA-2 RPKM', ebna2_rpkm),
    ('LMP-1 RPKM', lmp1_rpkm),
    ('log2(Lytic IE / Latent mRNA)', ie_lat_rat),
    ('WGS log2(CN) in MAGE (n=717)', lcn_m)
]:
    r_raw = pearson(x_ir1_raw, vec)
    vec_wb = demean(demean(vec, pops_m), bats_m)
    r_wb  = pearson(x_ir1_wb, vec_wb)
    print(f'  {name:<40}: marginal r = {r_raw:+.4f} (p={pv(r_raw,len(m_ir1)):.2e}) | within-pop+batch r = {r_wb:+.4f} (p={pv(r_wb,len(m_ir1)):.2e})')

# Also check IR1 vs log2(CN) across ALL 2,449 clean WGS lines (marginal vs within-pop)
w_ir1 = [s for s in clean_wgs if s in ir1]
x_w = [ir1[s] for s in w_ir1]; y_w = [math.log2(cn[s]) for s in w_ir1]
pops_w = [ped[s][0] for s in w_ir1]
r_w_raw = pearson(x_w, y_w)
r_w_wp  = pearson(demean(x_w, pops_w), demean(y_w, pops_w))
print(f'  WGS log2(CN) in ALL 1kGP (n={len(w_ir1)})      : marginal r = {r_w_raw:+.4f} (p={pv(r_w_raw,len(w_ir1)):.2e}) | within-26-pops   r = {r_w_wp:+.4f} (p={pv(r_w_wp,len(w_ir1)):.2e})')

# Check ONT IR1 replication of the cis-eQTL in MAGE!
m_ont_ir1 = [d for d in m_ir1 if d in ont_ir1]
xo_raw = [ont_ir1[d] for d in m_ont_ir1]
xo_wb  = demean(demean(xo_raw, [ped[d][0] for d in m_ont_ir1]), [by_d[d][0]['batch'] for d in m_ont_ir1])
yo_w1  = [math.log2((d_metric(d, lambda r: float(r['ir1_mapq0_reads'])*1e6/float(r['total_reads'])) + 0.1) / (d_metric(d, lambda r: float(r['latent_mrna_rpkm'])) + 0.1)) for d in m_ont_ir1]
yo_wb  = demean(demean(yo_w1, [ped[d][0] for d in m_ont_ir1]), [by_d[d][0]['batch'] for d in m_ont_ir1])
print(f'  [ONT IR1 validation, n={len(m_ont_ir1)}] ONT_IR1 vs log2(IR1_RNA / Latent_mRNA): marginal r = {pearson(xo_raw, yo_w1):+.4f} | within-pop+batch r = {pearson(xo_wb, yo_wb):+.4f} (p={pv(pearson(xo_wb, yo_wb), len(m_ont_ir1)):.2e})')

# ==============================================================================
# DESIGN B: Within-Population Viral CN-QTLs & Viral Cis/Trans eQTLs
# ==============================================================================
# Load LoFreq calls across all 2,502 samples
calls_by_site = defaultdict(dict) # (pos, ref, alt) -> {sample: vaf}
for line in open(f'{BASE}/variants/all_calls.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5:
        s, pos, ref, alt, vaf = f[0], int(f[1]), f[2], f[3], float(f[4])
        if s in clean_wgs:
            calls_by_site[(pos, ref, alt)][s] = vaf

# Exclude the 14 systematic artifact sites (>50% frequency or in 7638-8028 / 36416)
poly_sites = []
for (pos, ref, alt), sdict in calls_by_site.items():
    if 7600 <= pos <= 8035 or pos == 36416: continue
    n_clean = len(sdict)
    n_mage  = sum(1 for d in clean_mage if d in sdict)
    if 8 <= n_mage <= 600 and n_clean <= 1200:
        poly_sites.append((pos, ref, alt, n_clean, n_mage, sdict))

print(f'\n=== AIM 2.2 (DESIGN B): WITHIN-POPULATION VIRAL CN-QTLs & VIRAL CIS/TRANS eQTLs ({len(poly_sites)} polymorphic non-artifact sites with >=8 MAGE carriers) ===')
# Map each site to our 84 curated EBV genes
gene_intervals = [
    ('LMP-2B_ex','LATENT_MRNA',59,1680),('BNRF1','LYTIC_LATE',1736,5692),('EBER-1','POL3_EBER',6629,6795),('EBER-2','POL3_EBER',6956,7128),
    ('BCRF1','LYTIC_LATE',9675,10187),('EBNA-2','LATENT_MRNA',36216,37679),('BHLF1','LYTIC_EARLY',38287,38380),('BHRF1','LYTIC_EARLY',42088,42663),
    ('BFLF2','LYTIC_EARLY',43691,44647),('BFLF1','LYTIC_EARLY',44660,46235),('BFRF1A','LYTIC_EARLY',46236,46602),('BFRF1','LYTIC_EARLY',46603,47519),
    ('BFRF2','LYTIC_EARLY',47520,49218),('BFRF3','LYTIC_LATE',49219,49749),('BPLF1','LYTIC_LATE',49790,59231),('BOLF1','LYTIC_LATE',59240,62949),
    ('BORF1','LYTIC_LATE',62950,64044),('BORF2','LYTIC_EARLY',64119,66599),('BaRF1','LYTIC_EARLY',66612,67520),('BMRF1','LYTIC_EARLY',67611,68825),
    ('BMRF2','LYTIC_EARLY',68830,69903),('BSLF2_BMLF1','LYTIC_EARLY',70400,71940),('BSLF1','LYTIC_EARLY',71969,74593),('BSRF1','LYTIC_LATE',74636,75292),
    ('BLLF3','LYTIC_EARLY',75350,76186),('BLRF1','LYTIC_LATE',76259,76567),('BLRF2','LYTIC_LATE',76637,77125),('BLLF1','LYTIC_LATE',77142,79865),
    ('EBNA-3A','LATENT_MRNA',79955,82877),('EBNA-3B_3C','LATENT_MRNA',83065,89135),('BZLF2','LYTIC_LATE',89157,89828),('BZLF1','LYTIC_IE',89838,90943),
    ('BRLF1','LYTIC_IE',91078,92894),('BRRF1','LYTIC_EARLY',92895,93826),('BRRF2','LYTIC_LATE',94014,95627),('EBNA-1_N','LATENT_MRNA',95662,95950),
    ('EBNA-1_C','LATENT_MRNA',96600,97587),('BKRF2','LYTIC_LATE',97670,98064),('BKRF3','LYTIC_EARLY',98084,98832),('BKRF4','LYTIC_LATE',98846,99499),
    ('BBLF4','LYTIC_EARLY',99542,101915),('BBRF1','LYTIC_LATE',101916,103659),('BBRF2','LYTIC_LATE',103660,104496),('BBLF2_BBLF3','LYTIC_EARLY',105098,106800),
    ('BBRF3','LYTIC_LATE',106849,108066),('BBLF1','LYTIC_LATE',108459,108640),('BGLF5','LYTIC_EARLY',108687,110036),('BGLF4','LYTIC_EARLY',110037,111204),
    ('BGLF3.5','LYTIC_EARLY',111205,111652),('BGLF3','LYTIC_EARLY',111653,112649),('BGRF1_BDRF1','LYTIC_LATE',112650,113574),('BGLF2','LYTIC_LATE',113575,114562),
    ('BGLF1','LYTIC_LATE',114563,116055),('BDLF4','LYTIC_EARLY',116087,116711),('BDLF3.5','LYTIC_EARLY',116712,116945),('BDLF3','LYTIC_LATE',118074,118778),
    ('BDLF2','LYTIC_LATE',118839,120101),('BDLF1','LYTIC_LATE',120112,121017),('BcLF1','LYTIC_LATE',121033,125178),('BcRF1','LYTIC_EARLY',125179,127415),
    ('BTRF1','LYTIC_LATE',127416,128626),('BXLF2','LYTIC_LATE',128627,130747),('BXLF1','LYTIC_EARLY',130749,132570),('BXRF1','LYTIC_LATE',132571,133126),
    ('BVRF1','LYTIC_LATE',133127,134808),('BVLF1','LYTIC_EARLY',134840,135627),('BVRF2_BdRF1','LYTIC_LATE',135638,137455),('BILF2','LYTIC_LATE',137490,138236),
    ('LF3','B958_DEL',140570,143344),('LF2','B958_DEL',149035,150284),('LF1','B958_DEL',150285,151400),('BART_RPMS1','LATENT_BART',151551,153200),
    ('BILF1','LYTIC_EARLY',151703,152641),('BALF5','LYTIC_EARLY',153241,156288),('BALF4','LYTIC_LATE',156291,158850),('BALF3','LYTIC_EARLY',158865,159120),
    ('BART_BARF0','LATENT_BART',159121,160536),('BALF2','LYTIC_EARLY',160926,164312),('BALF1','LYTIC_EARLY',164397,164945),('BARF1','LYTIC_EARLY',165046,165711),
    ('LMP-2A_ex1','LATENT_MRNA',166103,166458),('BNLF2b','LYTIC_EARLY',166540,166836),('BNLF2a','LYTIC_EARLY',166846,167028),('LMP-1','LATENT_MRNA',168163,169474)
]
def pos_to_gene(p):
    for g, k, s, e in gene_intervals:
        if s <= p <= e: return g, k
    return 'intergenic', 'NA'

pops_wgs_all = [ped[s][0] for s in clean_wgs]
lcn_wgs_all  = [math.log2(cn[s]) for s in clean_wgs]
lcn_wgs_wp   = demean(lcn_wgs_all, pops_wgs_all)

pops_hold = [ped[s][0] for s in non_mage_wgs]
lcn_hold  = [math.log2(cn[s]) for s in non_mage_wgs]
lcn_hold_wp = demean(lcn_hold, pops_hold)

pops_mage = [ped[d][0] for d in clean_mage]
bats_mage = [by_d[d][0]['batch'] for d in clean_mage]
lcn_mage  = [math.log2(cn[d]) for d in clean_mage]
lcn_mage_wb = demean(demean(lcn_mage, pops_mage), bats_mage)
ie_mage   = [math.log2(d_metric(d, lambda r: float(r['lytic_ie_rpkm'])) + 0.1) for d in clean_mage]
ie_mage_wb = demean(demean(ie_mage, pops_mage), bats_mage)

class_col = {'LATENT_MRNA': 'latent_mrna_rpkm', 'LATENT_BART': 'latent_bart_rpkm', 'POL3_EBER': 'pol3_eber_rpkm',
             'LYTIC_IE': 'lytic_ie_rpkm', 'LYTIC_EARLY': 'lytic_early_rpkm', 'LYTIC_LATE': 'lytic_late_rpkm', 'B958_DEL': 'b958_del_rpkm'}

print(f'{"Site":<18} {"Gene":<12} {"n_WGS":>5} {"n_MAGE":>6} {"r_CN_raw":>8} {"r_CN_wp(2449)":>14} {"r_CN_hold(1732)":>15} {"r_IE_wb(717)":>13} {"r_cis_wb(717)":>14}')
for pos, ref, alt, n_clean, n_mage, sdict in sorted(poly_sites, key=lambda x: x[3], reverse=True):
    gname, kclass = pos_to_gene(pos)
    c_wgs = [1.0 if s in sdict else 0.0 for s in clean_wgs]
    c_wgs_wp = demean(c_wgs, pops_wgs_all)
    r_cn_raw = pearson(c_wgs, lcn_wgs_all)
    r_cn_wp  = pearson(c_wgs_wp, lcn_wgs_wp)

    c_hold = [1.0 if s in sdict else 0.0 for s in non_mage_wgs]
    r_cn_hold = pearson(demean(c_hold, pops_hold), lcn_hold_wp) if sum(c_hold) >= 3 else 0.0

    c_mage = [1.0 if d in sdict else 0.0 for d in clean_mage]
    c_mage_wb = demean(demean(c_mage, pops_mage), bats_mage)
    r_ie_wb = pearson(c_mage_wb, ie_mage_wb)

    # Cis-gene relative expression: log2((gene_rpkm + 0.1) / (class_rpkm + 0.1))
    r_cis_wb = 0.0
    if gname in gene_class and kclass in class_col:
        cis_vec = [math.log2((donor_gene[d][gname] + 0.1) / (d_metric(d, lambda r: float(r[class_col[kclass]])) + 0.1)) for d in clean_mage]
        r_cis_wb = pearson(c_mage_wb, demean(demean(cis_vec, pops_mage), bats_mage))

    print(f'{f"{pos}:{ref}>{alt}":<18} {gname:<12} {n_clean:5d} {n_mage:6d} {r_cn_raw:+8.4f} {r_cn_wp:+14.4f} {r_cn_hold:+15.4f} {r_ie_wb:+13.4f} {r_cis_wb:+14.4f}')

# ==============================================================================
# DESIGN D: Rare Regulatory Variant Burden & Wild-Type LF2 Co-Infection
# ==============================================================================
print('\n=== AIM 2.4 (DESIGN D1): RARE VIRAL PROMOTER/REPRESSOR VARIANT BURDEN (<=5 carriers, outside artifact mask) ===')
reg_windows = [
    ('Zp (BZLF1 promoter, 90944-91500)', 90944, 91500),
    ('Rp (BRLF1 promoter, 92896-93500)', 92896, 93500),
    ('oriP (FR+DS, 7315-7600 & 8050-9312)', 8050, 9312),
    ('EBNA-1 DBD (96600-97587)', 96600, 97587),
    ('oriLyt (40300-42000 & 52500-54000)', 40300, 42000),
    ('LMP-1 promoter/coding (168163-169600)', 168163, 169600),
]
lcn_wp_map = dict(zip(clean_wgs, lcn_wgs_wp))
ie_wb_map  = dict(zip(clean_mage, ie_mage_wb))
for wname, ws, we in reg_windows:
    carriers = set()
    for (pos, ref, alt), sdict in calls_by_site.items():
        if ws <= pos <= we and 1 <= len(sdict) <= 5:
            carriers.update(sdict.keys())
    c_wgs_list = [lcn_wp_map[s] for s in carriers if s in lcn_wp_map]
    c_mage_list = [ie_wb_map[s] for s in carriers if s in ie_wb_map]
    print(f'  {wname:<38}: n_WGS_carriers={len(c_wgs_list):3d}, mean within-pop log2(CN) shift = {(st.mean(c_wgs_list) if c_wgs_list else 0.0):+.3f} | n_MAGE={len(c_mage_list):2d}, mean IE shift = {(st.mean(c_mage_list) if c_mage_list else 0.0):+.3f}')

print(f'\n=== AIM 2.4 (DESIGN D2): WILD-TYPE EBV CO-INFECTION & LF2 REPRESSION TEST (n={len(coinf_mage)} co-infected vs {len(clean_mage)} B95-8 MAGE donors) ===')
lf2_clean = [donor_gene[d]['LF2'] for d in clean_mage]
lf2_coinf = [donor_gene[d]['LF2'] for d in coinf_mage]
lf1_clean = [donor_gene[d]['LF1'] for d in clean_mage]
lf1_coinf = [donor_gene[d]['LF1'] for d in coinf_mage]
# Late / IE progression ratio: log2((lytic_late_rpkm + 0.1) / (lytic_ie_rpkm + 0.1))
prog_clean = [math.log2((d_metric(d, lambda r: float(r['lytic_late_rpkm'])) + 0.1) / (d_metric(d, lambda r: float(r['lytic_ie_rpkm'])) + 0.1)) for d in clean_mage]
prog_coinf = [math.log2((d_metric(d, lambda r: float(r['lytic_late_rpkm'])) + 0.1) / (d_metric(d, lambda r: float(r['lytic_ie_rpkm'])) + 0.1)) for d in coinf_mage]
ie_clean_v = [d_metric(d, lambda r: float(r['lytic_ie_over_latent_log2'])) for d in clean_mage]
ie_coinf_v = [d_metric(d, lambda r: float(r['lytic_ie_over_latent_log2'])) for d in coinf_mage]
print(f'  LF1 RPKM (unique B95-8 del): clean median = {st.median(lf1_clean):.3f} (mean {st.mean(lf1_clean):.3f}) vs co-infected median = {st.median(lf1_coinf):.3f} (mean {st.mean(lf1_coinf):.3f})')
print(f'  LF2 RPKM (Rta inhibitor)   : clean median = {st.median(lf2_clean):.3f} (mean {st.mean(lf2_clean):.3f}) vs co-infected median = {st.median(lf2_coinf):.3f} (mean {st.mean(lf2_coinf):.3f})')
print(f'  log2(Lytic IE / Latent)    : clean median = {st.median(ie_clean_v):+.3f} (mean {st.mean(ie_clean_v):+.3f}) vs co-infected median = {st.median(ie_coinf_v):+.3f} (mean {st.mean(ie_coinf_v):+.3f})')
print(f'  log2(Lytic Late / Lytic IE): clean median = {st.median(prog_clean):+.3f} (mean {st.mean(prog_clean):+.3f}) vs co-infected median = {st.median(prog_coinf):+.3f} (mean {st.mean(prog_coinf):+.3f})')

with open(f'{BASE}/metadata/viral_qtl_scan_112.tsv', 'w') as out:
    out.write('site\tpos\tref\talt\tgene\tkinetic_class\tn_wgs_2449\tn_mage_717\tr_cn_raw\tr_cn_within_pop_2449\tr_cn_holdout_1732\tr_ie_within_pop_batch_717\tr_cis_within_pop_batch_717\n')
    for pos, ref, alt, n_clean, n_mage, sdict in sorted(poly_sites, key=lambda x: x[3], reverse=True):
        gname, kclass = pos_to_gene(pos)
        c_wgs = [1.0 if s in sdict else 0.0 for s in clean_wgs]
        r_cn_raw = pearson(c_wgs, lcn_wgs_all)
        r_cn_wp  = pearson(demean(c_wgs, pops_wgs_all), lcn_wgs_wp)
        c_hold = [1.0 if s in sdict else 0.0 for s in non_mage_wgs]
        r_cn_hold = pearson(demean(c_hold, pops_hold), lcn_hold_wp) if sum(c_hold) >= 3 else 0.0
        c_mage = [1.0 if d in sdict else 0.0 for d in clean_mage]
        c_mage_wb = demean(demean(c_mage, pops_mage), bats_mage)
        r_ie_wb = pearson(c_mage_wb, ie_mage_wb)
        r_cis_wb = 0.0
        if gname in gene_class and kclass in class_col:
            cis_vec = [math.log2((donor_gene[d][gname] + 0.1) / (d_metric(d, lambda r: float(r[class_col[kclass]])) + 0.1)) for d in clean_mage]
            r_cis_wb = pearson(c_mage_wb, demean(demean(cis_vec, pops_mage), bats_mage))
        out.write(f'{pos}:{ref}>{alt}\t{pos}\t{ref}\t{alt}\t{gname}\t{kclass}\t{n_clean}\t{n_mage}\t{r_cn_raw:.4f}\t{r_cn_wp:.4f}\t{r_cn_hold:.4f}\t{r_ie_wb:.4f}\t{r_cis_wb:.4f}\n')
