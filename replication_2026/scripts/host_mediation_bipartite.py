#!/usr/bin/env python3
"""Aim 3: Viral Gene Expression as a Causal Mediator of Host Gene Expression (n = 717 clean MAGE donors x 20,154 host genes)."""
import csv, gzip, math, statistics as st
from collections import defaultdict
import numpy as np

BASE = '/work/nfp8/EBV'
RU = 22737.0 / 3072.0

# 1. Load ENSG -> Gene Symbol map from human_nr_tx.fa
ensg2sym = {}
for line in open(f'{BASE}/metadata/human_nr_tx.fa'):
    if line.startswith('>'):
        p = line[1:].strip().split('|')
        if len(p) >= 3:
            ensg_base = p[1].split('.')[0]
            ensg2sym[p[1]] = p[2]
            ensg2sym[ensg_base] = p[2]

ped = {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    f = line.split()
    if len(f) >= 7 and f[1] != 'SampleID':
        ped[f[1]] = (int(f[4]), f[5], f[6]) # sex, pop, superpop

cn, wt, ir1, gcb = {}, {}, {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5:
        try: ir1[f[0]] = float(f[4]) * RU
        except ValueError: pass
for line in open(f'{BASE}/metadata/all_gcbias_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5 and f[4] != 'NA':
        try: gcb[f[0]] = float(f[4])
        except ValueError: pass

rna_rows = [r for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_rna_779.tsv'), delimiter='\t') if r['donor'] in cn]
by_d = defaultdict(list)
for r in rna_rows: by_d[r['donor']].append(r)

run_gene = defaultdict(dict)
for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_genes_779.tsv'), delimiter='\t'):
    run_gene[r['run']][r['gene']] = float(r['rpkm'])

def d_metric(d, fn):
    return st.mean([fn(r) for r in by_d[d]])
def d_gene(d, g):
    return st.mean([run_gene[r['run']][g] for r in by_d[d]])

# Load MAGE host expression matrix (20,154 genes x 731 donors)
BED = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/expression_quants/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz'
with gzip.open(BED, 'rt') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    s_all = hdr[4:]
    clean_idx = [(i, s) for i, s in enumerate(s_all) if s in by_d and wt.get(s, 0) <= 0.01 and s in ped and s in ir1]
    c_cols = [i for i, _ in clean_idx]
    donors = [s for _, s in clean_idx]
    gene_ids, gene_syms, Y_list = [], [], []
    for line in fh:
        f = line.rstrip('\n').split('\t')
        gid = f[3]
        gene_ids.append(gid)
        gene_syms.append(ensg2sym.get(gid, ensg2sym.get(gid.split('.')[0], gid)))
        Y_list.append([float(f[4+i]) for i in c_cols])

Y = np.array(Y_list, dtype=np.float64) # (G, N)
G, N = Y.shape

# Build covariate design matrix C: intercept + sex + 17 batch dummies + 26 pop dummies
sex_v = [float(ped[d][0]) for d in donors]
bats  = [by_d[d][0]['batch'] for d in donors]
pops  = [ped[d][1] for d in donors]
u_bats = sorted(set(bats))[1:]
u_pops = sorted(set(pops))[1:]

C_cols = [[1.0]*N, sex_v]
for b in u_bats: C_cols.append([1.0 if x == b else 0.0 for x in bats])
for p in u_pops: C_cols.append([1.0 if x == p else 0.0 for x in pops])
C = np.array(C_cols, dtype=np.float64).T # (N, K)
df_resid = N - C.shape[1] - 1

# Project out covariates C from Y and from any predictor vector x: M = I - C (C^T C)^+ C^T
P_C = C @ np.linalg.pinv(C)
def residualize_mat(M): # M is (G, N)
    return M - (M @ P_C)
def residualize_vec(v): # v is (N,)
    v = np.array(v, dtype=np.float64)
    v_r = v - (P_C @ v)
    return v_r / np.linalg.norm(v_r)

Y_res = residualize_mat(Y)
Y_norms = np.linalg.norm(Y_res, axis=1, keepdims=True)
Y_unit = Y_res / Y_norms

# Predictor vectors across the N = 717 clean donors
v_cn     = [math.log2(cn[d]) for d in donors]
v_ie     = [math.log2(d_metric(d, lambda r: float(r['lytic_ie_rpkm'])) + 0.1) for d in donors]
v_bz_spl = [math.log2(d_metric(d, lambda r: float(r['bzlf1_spliced_reads'])*1e6/float(r['total_reads'])) + 0.05) for d in donors]
v_bglf5  = [math.log2(d_gene(d, 'BGLF5') + 0.1) for d in donors]
v_late   = [math.log2(d_metric(d, lambda r: float(r['lytic_late_rpkm'])) + 0.1) for d in donors]
v_ie_lat = [d_metric(d, lambda r: float(r['lytic_ie_over_latent_log2'])) for d in donors]
v_lt_lat = [d_metric(d, lambda r: float(r['lytic_late_over_latent_log2'])) for d in donors]
v_lat    = [math.log2(d_metric(d, lambda r: float(r['latent_mrna_rpkm'])) + 0.1) for d in donors]
v_lmp1   = [math.log2(d_gene(d, 'LMP-1') + 0.1) for d in donors]
v_ebna2  = [math.log2(d_gene(d, 'EBNA-2') + 0.1) for d in donors]
v_bart   = [math.log2(d_metric(d, lambda r: float(r['latent_bart_rpkm'])) + 0.1) for d in donors]
v_eber   = [math.log2(d_metric(d, lambda r: float(r['pol3_eber_rpkm'])) + 0.01) for d in donors]
v_ir1    = [ir1[d] for d in donors]

def calc_stats(vec):
    u = residualize_vec(vec)
    r = Y_unit @ u # partial correlations across all 20,154 genes!
    t = np.abs(r) * np.sqrt(df_resid / np.maximum(1e-12, 1.0 - r*r))
    # Normal/Student-t p-value approximation (df=672)
    # Exact 2-sided Student-t p-value via incomplete beta or high-df normal+cornish-fisher
    z = t * (1.0 - 1.0 / (4.0 * df_resid)) / np.sqrt(1.0 + t*t / (2.0 * df_resid))
    erfc = np.vectorize(math.erfc)
    p = erfc(z / np.sqrt(2.0))
    # BH-FDR
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, G + 1)
    q = p * G / ranks
    q_min = np.minimum.accumulate(q[order[::-1]])[::-1]
    q_final = np.empty_like(q)
    q_final[order] = q_min
    pi0 = min(1.0, np.mean(p > 0.5) / 0.5)
    return r, p, q_final, 1.0 - pi0

print(f'=== AIM 3.1: SAME-VIAL VIRAL RNA vs 2019 WGS CN ACROSS {G} HUMAN GENES (n={N} clean donors, controlling sex+17 batches+26 pops, df={df_resid}) ===')
print(f'{"Viral Predictor":<34} {"Storey_pi1":>10} {"q<0.05_genes":>14} {"q<0.01_genes":>14} {"Mean_|r|":>9} {"Align_w_CN_r":>13}')
r_cn, p_cn, q_cn, pi1_cn = calc_stats(v_cn)
stats_map = {}
for name, vec in [
    ('1. 2019 WGS log2(CN)', v_cn),
    ('2. Same-vial Lytic IE RPKM', v_ie),
    ('3. Spliced BZLF1 (Zta) RPM', v_bz_spl),
    ('4. BGLF5 (Host-Shutoff Nuclease)', v_bglf5),
    ('5. Lytic Late Structural RPKM', v_late),
    ('6. log2(Lytic IE / Latent mRNA)', v_ie_lat),
    ('7. log2(Lytic Late / Latent mRNA)', v_lt_lat),
    ('8. Latent Pol II mRNA RPKM', v_lat),
    ('9. LMP-1 (NF-kB Latent Driver)', v_lmp1),
    ('10. EBNA-2 (MYC Latent Driver)', v_ebna2),
    ('11. BART lncRNA RPKM', v_bart),
    ('12. Pol III EBER1/2 RPKM', v_eber),
    ('13. Exogenous IR1 Repeat Count', v_ir1),
]:
    r_v, p_v, q_v, pi1_v = calc_stats(vec)
    stats_map[name] = (r_v, p_v, q_v)
    n05 = int(np.sum(q_v < 0.05))
    n01 = int(np.sum(q_v < 0.01))
    align_cn = float(np.corrcoef(r_cn, r_v)[0, 1])
    print(f'{name:<34} {pi1_v*100:9.1f}% {f"{n05} ({n05/G*100:.1f}%)":>14} {f"{n01} ({n01/G*100:.1f}%)":>14} {float(np.mean(np.abs(r_v))):9.4f} {align_cn:+13.4f}')

# Formal Baron-Kenny / Sobel Mediation: Conditioning WGS log2(CN) on Same-Vial Viral Lytic RNA (BGLF5 + IE + Late + BART)
# And conditioning Viral Lytic RNA on 2019 WGS log2(CN)!
u_cn    = residualize_vec(v_cn)
u_bglf5 = residualize_vec(v_bglf5)
u_ie    = residualize_vec(v_ie)
u_late  = residualize_vec(v_late)

# Partial correlation of WGS_CN with host genes controlling for BGLF5 + IE + Late
L_mat = np.column_stack([C, v_bglf5, v_ie, v_late, v_bz_spl])
P_L   = L_mat @ np.linalg.pinv(L_mat)
Y_res_L = Y - (Y @ P_L)
Y_unit_L = Y_res_L / np.linalg.norm(Y_res_L, axis=1, keepdims=True)
u_cn_L = np.array(v_cn) - (P_L @ np.array(v_cn))
u_cn_L = u_cn_L / np.linalg.norm(u_cn_L)
r_cn_given_lytic = Y_unit_L @ u_cn_L
t_cn_L = np.abs(r_cn_given_lytic) * np.sqrt((df_resid - 4) / np.maximum(1e-12, 1.0 - r_cn_given_lytic**2))
erfc = np.vectorize(math.erfc)
p_cn_L = erfc(t_cn_L / np.sqrt(2.0))
order_L = np.argsort(p_cn_L); ranks_L = np.empty_like(order_L); ranks_L[order_L] = np.arange(1, G + 1)
q_cn_L = np.empty_like(p_cn_L); q_cn_L[order_L] = np.minimum.accumulate((p_cn_L * G / ranks_L)[order_L[::-1]])[::-1]

# Conversely, partial correlation of BGLF5 with host genes controlling for 2019 WGS_CN!
C_cn = np.column_stack([C, v_cn])
P_Ccn = C_cn @ np.linalg.pinv(C_cn)
Y_res_cn = Y - (Y @ P_Ccn)
Y_unit_cn = Y_res_cn / np.linalg.norm(Y_res_cn, axis=1, keepdims=True)
u_bg_cn = np.array(v_bglf5) - (P_Ccn @ np.array(v_bglf5))
u_bg_cn = u_bg_cn / np.linalg.norm(u_bg_cn)
r_bglf5_given_cn = Y_unit_cn @ u_bg_cn
t_bg_cn = np.abs(r_bglf5_given_cn) * np.sqrt((df_resid - 1) / np.maximum(1e-12, 1.0 - r_bglf5_given_cn**2))
p_bg_cn = erfc(t_bg_cn / np.sqrt(2.0))
order_B = np.argsort(p_bg_cn); ranks_B = np.empty_like(order_B); ranks_B[order_B] = np.arange(1, G + 1)
q_bg_cn = np.empty_like(p_bg_cn); q_bg_cn[order_B] = np.minimum.accumulate((p_bg_cn * G / ranks_B)[order_B[::-1]])[::-1]

sig_cn = q_cn < 0.05
med_frac = 1.0 - np.mean(r_cn_given_lytic[sig_cn]**2) / np.mean(r_cn[sig_cn]**2)
print(f'\n=== AIM 3.1B: CAUSAL MEDIATION ASYMMETRY (2019 WGS CN vs SAME-VIAL 2023 LYTIC RNA) ===')
print(f'  WGS log2(CN) unadjusted for Lytic RNA : {int(np.sum(q_cn < 0.05)):5d} genes at q < 0.05 (mean R^2 on sig genes = {float(np.mean(r_cn[sig_cn]**2))*100:.2f}%)')
print(f'  WGS log2(CN) controlling for Lytic RNA: {int(np.sum(q_cn_L < 0.05)):5d} genes at q < 0.05 (mean R^2 on sig genes = {float(np.mean(r_cn_given_lytic[sig_cn]**2))*100:.2f}% -> {med_frac*100:.1f}% of variance MEDIATED by same-vial Lytic RNA!)')
print(f'  Same-vial BGLF5 controlling for WGS CN: {int(np.sum(q_bg_cn < 0.05)):5d} genes STILL significant at q < 0.05 even after removing 2019 WGS CN!')

# Aim 3.2: Inspect canonical host biological modules across the 4 orthogonal viral axes!
r_bglf5 = stats_map['4. BGLF5 (Host-Shutoff Nuclease)'][0]
r_ie_v  = stats_map['2. Same-vial Lytic IE RPKM'][0]
r_lmp1  = stats_map['9. LMP-1 (NF-kB Latent Driver)'][0]
r_ebna2 = stats_map['10. EBNA-2 (MYC Latent Driver)'][0]
r_bart  = stats_map['11. BART lncRNA RPKM'][0]
r_eber  = stats_map['12. Pol III EBER1/2 RPKM'][0]
r_ir1_v = stats_map['13. Exogenous IR1 Repeat Count'][0]

sym_idx = {}
for i, s in enumerate(gene_syms):
    if s not in sym_idx: sym_idx[s] = i

print('\n=== AIM 3.2: BIPARTITE VIRAL-TO-HOST FUNCTIONAL MODULE DISSECTION (Within-Pop+Batch Partial r across n=717) ===')
print(f'{"Module / Host Gene":<26} {"r(WGS_CN)":>10} {"r(Lytic_IE)":>12} {"r(BGLF5)":>10} {"r(LMP-1)":>10} {"r(EBNA-2)":>10} {"r(BART)":>9} {"r(IR1_DNA)":>11}')
marker_panels = [
    ('--- Plasmablast / ER / UPR ---', []),
    ('XBP1 (UPR / Zp activator)', 'XBP1'),
    ('PRDM1 (Blimp-1)', 'PRDM1'),
    ('IRF4 (Plasmablast TF)', 'IRF4'),
    ('MZB1 (ER chaperone)', 'MZB1'),
    ('JCHAIN (Ig J-chain)', 'JCHAIN'),
    ('HSPA5 (BiP / GRP78)', 'HSPA5'),
    ('SSR4 (ER translocon)', 'SSR4'),
    ('--- B-Cell Identity ---', []),
    ('PAX5 (B-cell master TF)', 'PAX5'),
    ('EBF1 (B-cell factor 1)', 'EBF1'),
    ('SPIB (ETS B-cell TF)', 'SPIB'),
    ('--- Cell Cycle / p53 Arrest ---', []),
    ('CDKN1A (p21 cip1)', 'CDKN1A'),
    ('GADD45A (Stress arrest)', 'GADD45A'),
    ('MKI67 (Proliferation)', 'MKI67'),
    ('CDK1 (Mitotic cyclin)', 'CDK1'),
    ('MYC (Proto-oncogene)', 'MYC'),
    ('--- NF-kB / Latent III Targets ---', []),
    ('TRAF1 (LMP1 target)', 'TRAF1'),
    ('BIRC3 (cIAP2 anti-apopt)', 'BIRC3'),
    ('NFKBIA (IkB-alpha)', 'NFKBIA'),
    ('CD40 (LMP1 co-receptor)', 'CD40'),
    ('ICAM1 (Adhesion)', 'ICAM1'),
    ('--- Type I IFN / Innate ISGs ---', []),
    ('ISG15 (Ubiquitin-like ISG)', 'ISG15'),
    ('MX1 (Antiviral GTPase)', 'MX1'),
    ('IFIT1 (IFN-induced)', 'IFIT1'),
    ('OAS1 (2-5A synthetase)', 'OAS1'),
    ('STAT1 (IFN TF)', 'STAT1'),
]
for label, sym in marker_panels:
    if not sym:
        print(f'  {label}')
        continue
    if sym in sym_idx:
        idx = sym_idx[sym]
        print(f'  {label:<24} {r_cn[idx]:+10.4f} {r_ie_v[idx]:+12.4f} {r_bglf5[idx]:+10.4f} {r_lmp1[idx]:+10.4f} {r_ebna2[idx]:+10.4f} {r_bart[idx]:+9.4f} {r_ir1_v[idx]:+11.4f}')

# Aim 3.3: Multi-Parameter Out-of-Sample Vector Alignment Test (24 non-IBS/GBR pops -> IBS vs GBR DE)
ibs_mask = np.array([ped[d][1] == 'IBS' for d in donors])
gbr_mask = np.array([ped[d][1] == 'GBR' for d in donors])
other_mask = ~(ibs_mask | gbr_mask)

# Sex+batch adjusted expression for IBS vs GBR mean difference
C_sb = np.array(C_cols[:2+len(u_bats)], dtype=np.float64).T
P_sb = C_sb @ np.linalg.pinv(C_sb)
Y_sb = Y - (Y @ P_sb)
de_ibs_gbr = np.mean(Y_sb[:, ibs_mask], axis=1) - np.mean(Y_sb[:, gbr_mask], axis=1)

# Fit within-population slopes strictly inside the 24 non-IBS/non-GBR populations!
C_24 = C[other_mask, :]
P_24 = C_24 @ np.linalg.pinv(C_24)
Y_24 = Y[:, other_mask] - (Y[:, other_mask] @ P_24)

def pred_24(preds):
    X_24 = np.column_stack([np.array(p)[other_mask] for p in preds])
    X_24_res = X_24 - (P_24 @ X_24)
    # Betas across all 20,154 genes estimated strictly in the 24 non-IBS/GBR populations:
    betas = Y_24 @ X_24_res @ np.linalg.pinv(X_24_res.T @ X_24_res) # (G, P)
    # Multiply by the observed viral difference between IBS and GBR:
    X_all_sb = np.column_stack([np.array(p) - (P_sb @ np.array(p)) for p in preds])
    delta_viral = np.mean(X_all_sb[ibs_mask, :], axis=0) - np.mean(X_all_sb[gbr_mask, :], axis=0)
    pred_de = betas @ delta_viral
    r_all = float(np.corrcoef(de_ibs_gbr, pred_de)[0, 1])
    r_sig = float(np.corrcoef(de_ibs_gbr[sig_cn], pred_de[sig_cn])[0, 1])
    return r_all, r_sig

print('\n=== AIM 3.3: MULTI-PARAMETER OUT-OF-SAMPLE VECTOR ALIGNMENT (Trained in 24 non-IBS/GBR pops -> Predicting IBS vs GBR DE) ===')
for name, p_list in [
    ('1. 2019 WGS log2(CN) alone', [v_cn]),
    ('2. Same-vial Lytic RNA (BGLF5 + IE + Late) alone', [v_bglf5, v_ie, v_late]),
    ('3. Combined Viral DNA CN + 4-Axis Viral RNA (CN + BGLF5 + IE + LMP1 + BART)', [v_cn, v_bglf5, v_ie, v_lmp1, v_bart]),
]:
    ra, rs = pred_24(p_list)
    print(f'  {name:<72}: All 20,154 genes r = {ra:+.4f} (R^2 = {ra**2*100:5.1f}%) | q<0.05 genes r = {rs:+.4f} (R^2 = {rs**2*100:5.1f}%)')

# Also compute 1-df conditioning on Spliced BZLF1 alone and export full 20,154-gene table
C_bz = np.column_stack([C, v_bz_spl])
P_Cbz = C_bz @ np.linalg.pinv(C_bz)
Y_res_bz = Y - (Y @ P_Cbz)
Y_unit_bz = Y_res_bz / np.linalg.norm(Y_res_bz, axis=1, keepdims=True)
u_cn_bz = np.array(v_cn) - (P_Cbz @ np.array(v_cn))
u_cn_bz = u_cn_bz / np.linalg.norm(u_cn_bz)
r_cn_given_bz = Y_unit_bz @ u_cn_bz

u_bz_cn = np.array(v_bz_spl) - (P_Ccn @ np.array(v_bz_spl))
u_bz_cn = u_bz_cn / np.linalg.norm(u_bz_cn)
r_bz_given_cn = Y_unit_cn @ u_bz_cn

X_24_cn = np.array(v_cn)[other_mask] - (P_24 @ np.array(v_cn)[other_mask])
beta_cn_24 = (Y_24 @ X_24_cn) / float(X_24_cn @ X_24_cn)

r_bz_v, p_bz_v, q_bz_v = stats_map['3. Spliced BZLF1 (Zta) RPM']
r_bg_v, p_bg_v, q_bg_v = stats_map['4. BGLF5 (Host-Shutoff Nuclease)']
r_lt_v, p_lt_v, q_lt_v = stats_map['7. log2(Lytic Late / Latent mRNA)']
r_lm_v, p_lm_v, q_lm_v = stats_map['9. LMP-1 (NF-kB Latent Driver)']
r_eb_v, p_eb_v, q_eb_v = stats_map['10. EBNA-2 (MYC Latent Driver)']
r_ba_v, p_ba_v, q_ba_v = stats_map['11. BART lncRNA RPKM']
r_ir_v, p_ir_v, q_ir_v = stats_map['13. Exogenous IR1 Repeat Count']

with open(f'{BASE}/metadata/host_viral_mediation_20154.tsv', 'w') as out:
    out.write('gene_id\tgene_symbol\tr_wgs_cn\tq_wgs_cn\tr_spliced_bzlf1\tq_spliced_bzlf1\tr_bglf5\tq_bglf5\tr_late_over_latent\tq_late_over_latent\tr_lmp1\tq_lmp1\tr_ebna2\tr_bart\tr_ir1_dna\tq_ir1_dna\tr_wgs_cn_given_bzlf1\tr_bzlf1_given_wgs_cn\tbeta_cn_24pops\tde_ibs_minus_gbr\n')
    for i in range(G):
        out.write(f'{gene_ids[i]}\t{gene_syms[i]}\t{r_cn[i]:.4f}\t{q_cn[i]:.3e}\t{r_bz_v[i]:.4f}\t{q_bz_v[i]:.3e}\t{r_bg_v[i]:.4f}\t{q_bg_v[i]:.3e}\t{r_lt_v[i]:.4f}\t{q_lt_v[i]:.3e}\t{r_lm_v[i]:.4f}\t{q_lm_v[i]:.3e}\t{r_eb_v[i]:.4f}\t{r_ba_v[i]:.4f}\t{r_ir_v[i]:.4f}\t{q_ir_v[i]:.3e}\t{r_cn_given_bz[i]:.4f}\t{r_bz_given_cn[i]:.4f}\t{beta_cn_24[i]:.4f}\t{de_ibs_gbr[i]:.4f}\n')
