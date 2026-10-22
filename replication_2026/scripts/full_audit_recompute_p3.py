import csv, gzip, math, os, random, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
RU = 22737 / 3072.0

pop, spop = {}, {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    f = line.split()
    if len(f) >= 7 and f[1] != 'SampleID':
        pop[f[1]] = f[5]; spop[f[1]] = f[6]

cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try:
            v = float(f[3]); w = float(f[9])
            if v > 0: cn[f[0]] = v; wt[f[0]] = w
        except ValueError: pass

unrel_2502 = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t') if r['sample'] in cn}
coinf = {s for s, v in wt.items() if v > 0.01}
clean_2502 = unrel_2502 - coinf

ont_cn = {}
for line in open(f'{BASE}/metadata/all_ont_904.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try:
            v = float(f[5])
            if v > 0: ont_cn[f[0]] = v
        except ValueError: pass

low_cn, low_ed = {}, {}
for line in open(f'{BASE}/metadata/all_lowcov_2499.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 6:
        try:
            ed = float(f[2]); v = float(f[5])
            low_ed[f[0]] = ed
            if v >= 0: low_cn[f[0]] = v
        except ValueError: pass

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    da = [x - ma for x in a]; db = [y - mb for y in b]
    den = math.sqrt(sum(x*x for x in da) * sum(y*y for y in db))
    return sum(x*y for x, y in zip(da, db)) / den if den else 0.0

# === 6. THREE-PLATFORM CN RELIABILITY AT FULL N ===
print('=== 6. THREE-PLATFORM CN COMPARISON AT FULL N ===')
n_blood_05 = sum(1 for v in low_ed.values() if v < 0.5)
n_gap = sum(1 for v in low_ed.values() if 0.46 <= v <= 19.65)
n_lcl_low = sum(1 for v in low_ed.values() if v >= 5.0)
print(f'  Low-coverage (n={len(low_ed)}): <0.5x (blood-derived) = {n_blood_05} ({100*n_blood_05/len(low_ed):.1f}%), in 0.46-19.65x gap = {n_gap}, >=5x (LCL) = {n_lcl_low}')

s_30_ont = sorted(s for s in clean_2502 if s in ont_cn and ont_cn[s] > 0.5)
s_30_low = sorted(s for s in clean_2502 if s in low_cn and low_ed[s] >= 5.0 and low_cn[s] > 0.5)
s_ont_low = sorted(s for s in s_30_ont if s in s_30_low)

r_30_ont = pearson([math.log2(cn[s]) for s in s_30_ont], [math.log2(ont_cn[s]) for s in s_30_ont])
r_30_low = pearson([math.log2(cn[s]) for s in s_30_low], [math.log2(low_cn[s]) for s in s_30_low])
r_ont_low = pearson([math.log2(ont_cn[s]) for s in s_ont_low], [math.log2(low_cn[s]) for s in s_ont_low])

print(f'  ONT vs 30x     (n={len(s_30_ont)} clean): r = {r_30_ont:.4f} (R^2 = {r_30_ont**2:.4f}), median ONT/30x = {st.median([ont_cn[s]/cn[s] for s in s_30_ont]):.4f}')
print(f'  Lowcov vs 30x  (n={len(s_30_low)} clean): r = {r_30_low:.4f} (R^2 = {r_30_low**2:.4f}), median lowcov/30x = {st.median([low_cn[s]/cn[s] for s in s_30_low]):.4f}')
print(f'  ONT vs Lowcov  (n={len(s_ont_low)} clean): r = {r_ont_low:.4f} (R^2 = {r_ont_low**2:.4f})')
# 3-indicator reliability on the 3-way overlap (and pairwise):
r12_3 = pearson([math.log2(cn[s]) for s in s_ont_low], [math.log2(ont_cn[s]) for s in s_ont_low])
r13_3 = pearson([math.log2(cn[s]) for s in s_ont_low], [math.log2(low_cn[s]) for s in s_ont_low])
r23_3 = r_ont_low
print(f'  3-indicator reliability (pairwise): 30x = {r_30_ont*r_30_low/r_ont_low:.4f}, ONT = {r_30_ont*r_ont_low/r_30_low:.4f}, Lowcov = {r_30_low*r_ont_low/r_30_ont:.4f}')
print(f'  3-indicator reliability (n={len(s_ont_low)} 3-way matched): 30x = {r12_3*r13_3/r23_3:.4f}, ONT = {r12_3*r23_3/r13_3:.4f}, Lowcov = {r13_3*r23_3/r12_3:.4f}')

# === 7. MAGE EXPRESSION ANOVA WITH CN ADDED + VECTOR ALIGNMENT TEST ===
print('\n=== 7. MAGE EXPRESSION ANOVA: HOW MUCH POPULATION VARIANCE IS ABSORBED BY EBV CN? ===')
import numpy as np

mage_meta = {}
for m in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sample.metadata.MAGE.v1.0.txt'), delimiter='\t'):
    mage_meta[m['sample_kgpID']] = m

ir1 = {}
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5 and f[4] != 'NA':
        try: ir1[f[0]] = float(f[4]) * RU
        except ValueError: pass

BED = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/expression_quants/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz'
with gzip.open(BED, 'rt') as fh:
    hdr = next(fh).rstrip('\n').split('\t')
    sids = hdr[4:]
    keep_idx = [j for j, s in enumerate(sids) if s in cn and s in mage_meta]
    keep_s   = [sids[j] for j in keep_idx]
    genes, Y_list = [], []
    for line in fh:
        p = line.rstrip('\n').split('\t')
        genes.append(p[3])
        Y_list.append([float(p[4 + j]) for j in keep_idx])
Y = np.array(Y_list, dtype=np.float64) # (G, N)
G, N = Y.shape
print(f'  Loaded expression matrix: {G} genes x {N} MAGE donors with CN')

# Build design matrices matching MAGE's generate_expression_ANOVAs.R:
# Model 0 (baseline): sex + batch
# Then sequential R^2 (Type I/II PVE) of continentalGroup and population BEFORE vs AFTER regressing out log2(CN)!
sex_v = [1.0 if mage_meta[s]['sex'] == 'F' else 0.0 for s in keep_s]
batches = sorted({mage_meta[s]['batch'] for s in keep_s})
pops    = sorted({mage_meta[s]['population'] for s in keep_s})
conts   = sorted({mage_meta[s]['continentalGroup'] for s in keep_s})
lcn_v   = np.array([math.log2(cn[s]) for s in keep_s], dtype=np.float64)

X_base = np.column_stack([[1.0]*N, sex_v] + [[1.0 if mage_meta[s]['batch']==b else 0.0 for s in keep_s] for b in batches[1:]])
X_cn   = np.column_stack([X_base, lcn_v])
D_pop  = np.column_stack([[1.0 if mage_meta[s]['population']==p else 0.0 for s in keep_s] for p in pops[1:]])
D_cont = np.column_stack([[1.0 if mage_meta[s]['continentalGroup']==c else 0.0 for s in keep_s] for c in conts[1:]])

def resid_matrix(Y_mat, X_mat):
    # Y_mat is (G, N), X_mat is (N, K)
    Q, _ = np.linalg.qr(X_mat)
    return Y_mat - (Y_mat @ Q) @ Q.T

def pve_from_resid(Y_res, D_mat, X_cov):
    # fraction of residual variance (or total variance) explained by D_mat after X_cov
    sst_res = np.sum(Y_res**2, axis=1)
    Y_after = resid_matrix(Y_res, np.column_stack([X_cov, D_mat]))
    ssr_after = np.sum(Y_after**2, axis=1)
    return (sst_res - ssr_after) / sst_res

Y_res_base = resid_matrix(Y, X_base)
Y_res_cn   = resid_matrix(Y, X_cn)

# Wait: how much variance of Y_res_base does log2(CN) explain by itself?
sst_base = np.sum(Y_res_base**2, axis=1)
sst_cn   = np.sum(Y_res_cn**2, axis=1)
pve_cn_alone = (sst_base - sst_cn) / sst_base

pve_pop_before  = pve_from_resid(Y_res_base, D_pop, X_base)
pve_pop_after   = (sst_cn - np.sum(resid_matrix(Y_res_cn, np.column_stack([X_cn, D_pop]))**2, axis=1)) / sst_base # as share of baseline residual variance
pve_cont_before = pve_from_resid(Y_res_base, D_cont, X_base)
pve_cont_after  = (sst_cn - np.sum(resid_matrix(Y_res_cn, np.column_stack([X_cn, D_cont]))**2, axis=1)) / sst_base

print(f'  Mean PVE across all {G} genes (as % of sex+batch residual variance):')
print(f'    EBV log2(CN) alone      : mean {100*np.mean(pve_cn_alone):.3f}%, median {100*np.median(pve_cn_alone):.3f}%, max {100*np.max(pve_cn_alone):.2f}%')
print(f'    Population (26) BEFORE  : mean {100*np.mean(pve_pop_before):.3f}% -> AFTER controlling CN: {100*np.mean(pve_pop_after):.3f}% (relative drop: {100*(1 - np.mean(pve_pop_after)/np.mean(pve_pop_before)):.1f}%)')
print(f'    Continental (5) BEFORE  : mean {100*np.mean(pve_cont_before):.3f}% -> AFTER controlling CN: {100*np.mean(pve_cont_after):.3f}% (relative drop: {100*(1 - np.mean(pve_cont_after)/np.mean(pve_cont_before)):.1f}%)')

# How about among the CN-associated genes (BH FDR q < 0.05 within-population)?
# Compute exact within-population OLS t-stat and p-value for log2(CN):
X_pop_cov = np.column_stack([X_base, D_pop])
Y_w = resid_matrix(Y, X_pop_cov)
x_w = resid_matrix(lcn_v.reshape(1, -1), X_pop_cov).ravel()
r_within = (Y_w @ x_w) / np.sqrt(np.sum(Y_w**2, axis=1) * np.sum(x_w**2))
dof_w = N - X_pop_cov.shape[1] - 1
t_within = r_within * np.sqrt(dof_w / np.maximum(1 - r_within**2, 1e-12))
from scipy.stats import t as t_dist
p_within = 2 * t_dist.sf(np.abs(t_within), dof_w)
# BH FDR
order = np.argsort(p_within)
q_within = np.empty_like(p_within)
prev = 1.0
for rk in range(G, 0, -1):
    idx = order[rk - 1]
    prev = min(prev, p_within[idx] * G / rk)
    q_within[idx] = prev
pi0_est = min(1.0, 2.0 * np.mean(p_within > 0.5))
pi1_est = 1.0 - pi0_est
sig_q05 = q_within < 0.05
sig_q01 = q_within < 0.01
print(f'\n  Rigorous Within-Population Association (controlling sex + batch + 26 populations, dof={dof_w}):')
print(f'    Storey pi_1 (true non-null fraction of transcriptome): {100*pi1_est:.1f}%')
print(f'    BH FDR q < 0.05 : {np.sum(sig_q05)} / {G} ({100*np.mean(sig_q05):.1f}%)')
print(f'    BH FDR q < 0.01 : {np.sum(sig_q01)} / {G} ({100*np.mean(sig_q01):.1f}%)')
print(f'    Among {np.sum(sig_q05)} CN-associated genes (q<0.05):')
print(f'      Pop PVE BEFORE: {100*np.mean(pve_pop_before[sig_q05]):.2f}% -> AFTER CN: {100*np.mean(pve_pop_after[sig_q05]):.2f}% (relative drop: {100*(1 - np.mean(pve_pop_after[sig_q05])/np.mean(pve_pop_before[sig_q05])):.1f}%)')
print(f'      Cont PVE BEFORE: {100*np.mean(pve_cont_before[sig_q05]):.2f}% -> AFTER CN: {100*np.mean(pve_cont_after[sig_q05]):.2f}% (relative drop: {100*(1 - np.mean(pve_cont_after[sig_q05])/np.mean(pve_cont_before[sig_q05])):.1f}%)')

# VECTOR ALIGNMENT TEST:
# Does a gene's WITHIN-population CN effect (beta_CN_within, totally orthogonal to population!)
# predict the BETWEEN-population expression difference between high-CN IBS and low-CN GBR (two European pops)?
ibs_mask = np.array([mage_meta[s]['population'] == 'IBS' for s in keep_s])
gbr_mask = np.array([mage_meta[s]['population'] == 'GBR' for s in keep_s])
diff_ibs_gbr = np.mean(Y_res_base[:, ibs_mask], axis=1) - np.mean(Y_res_base[:, gbr_mask], axis=1)
# Also compute beta_CN_within EXCLUDING IBS and GBR so there is ZERO sample overlap!
non_ig = ~(ibs_mask | gbr_mask)
X_pop_noig = np.column_stack([X_base[non_ig], D_pop[non_ig]])
Y_w_noig = resid_matrix(Y[:, non_ig], X_pop_noig)
x_w_noig = resid_matrix(lcn_v[non_ig].reshape(1, -1), X_pop_noig).ravel()
beta_cn_noig = (Y_w_noig @ x_w_noig) / np.sum(x_w_noig**2)

r_vec_all = np.corrcoef(beta_cn_noig, diff_ibs_gbr)[0, 1]
r_vec_sig = np.corrcoef(beta_cn_noig[sig_q05], diff_ibs_gbr[sig_q05])[0, 1]
print(f'\n  SMOKING-GUN VECTOR ALIGNMENT TEST (Out-of-sample within-pop beta_CN [24 non-IBS/GBR pops] vs [IBS - GBR] expression shift):')
print(f'    Across all {G} genes               : r = {r_vec_all:+.4f}')
print(f'    Across {np.sum(sig_q05)} CN-associated genes (q<0.05): r = {r_vec_sig:+.4f}')
