#!/usr/bin/env python3
"""Aim 1: High-resolution dissection of all 84 EBV genes across 779 MAGE libraries:
  1. Technical (n=24) vs Cross-thaw biological (n=46) repeatability per viral gene.
  2. Within-population+batch association with log2(CN) across 717 clean donors.
  3. Partial correlation with log2(CN) conditioning on spliced BZLF1 (abortive vs productive lytic replication).
  4. 3'-coterminal poly(A) cluster step-up flux analysis.
"""
import csv, math, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try: cn[f[0]] = float(f[3]); wt[f[0]] = float(f[9])
        except ValueError: pass

rna_rows = {r['run']: r for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_rna_779.tsv'), delimiter='\t') if r['donor'] in cn}

# Load 84 genes x 779 runs
run_gene = defaultdict(dict)
gene_meta = {}
for r in csv.DictReader(open(f'{BASE}/metadata/all_mage_ebv_genes_779.tsv'), delimiter='\t'):
    run = r['run']
    if run not in rna_rows: continue
    g = r['gene']
    run_gene[run][g] = float(r['rpkm'])
    if g not in gene_meta:
        gene_meta[g] = (r['kinetic_class'], int(r['length_bp']))

genes = list(gene_meta.keys())
by_d_runs = defaultdict(list)
for run, r in rna_rows.items():
    by_d_runs[r['donor']].append(run)

mult_d = {d: rs for d, rs in by_d_runs.items() if len(rs) >= 2}
same_pairs, diff_pairs = [], []
for d, rs in mult_d.items():
    for i in range(len(rs)):
        for j in range(i+1, len(rs)):
            if rna_rows[rs[i]]['batch'] == rna_rows[rs[j]]['batch']:
                same_pairs.append((rs[i], rs[j]))
            else:
                diff_pairs.append((rs[i], rs[j]))

def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx)**2 for a in x) * sum((b - my)**2 for b in y))
    return num / den if den > 0 else 0.0

def demean(vals, groups):
    gm = defaultdict(list)
    for v, g in zip(vals, groups): gm[g].append(v)
    m = {g: st.mean(vs) for g, vs in gm.items()}
    return [v - m[g] for v, g in zip(vals, groups)]

def resid(y, x):
    mx, my = st.mean(x), st.mean(y)
    vx = sum((a - mx)**2 for a in x)
    if vx <= 0: return y
    b = sum((a - mx)*(b - my) for a, b in zip(x, y)) / vx
    return [b_i - (my + b*(a_i - mx)) for a_i, b_i in zip(x, y)]

clean_donors = sorted([d for d in by_d_runs if wt.get(d, 0) <= 0.01])
pops = [rna_rows[by_d_runs[d][0]]['population'] for d in clean_donors]
bats = [rna_rows[by_d_runs[d][0]]['batch'] for d in clean_donors]
lcn = [math.log2(cn[d]) for d in clean_donors]
lcn_wb = demean(demean(lcn, pops), bats)

bz_spl = [math.log2(st.mean([float(rna_rows[r]['bzlf1_spliced_reads'])*1e6/float(rna_rows[r]['total_reads']) for r in by_d_runs[d]]) + 0.05) for d in clean_donors]
bz_wb  = demean(demean(bz_spl, pops), bats)
lcn_res_bz = resid(lcn_wb, bz_wb)

results = []
for g in genes:
    kclass, glen = gene_meta[g]
    # Replicate correlations (log2(rpkm + 0.1))
    r_tech = pearson([math.log2(run_gene[r1][g]+0.1) for r1, _ in same_pairs],
                     [math.log2(run_gene[r2][g]+0.1) for _, r2 in same_pairs])
    r_bio  = pearson([math.log2(run_gene[r1][g]+0.1) for r1, _ in diff_pairs],
                     [math.log2(run_gene[r2][g]+0.1) for _, r2 in diff_pairs])
    # Donor-level mean RPKM
    d_rpkm = [st.mean([run_gene[r][g] for r in by_d_runs[d]]) for d in clean_donors]
    d_log  = [math.log2(v + 0.1) for v in d_rpkm]
    d_wb   = demean(demean(d_log, pops), bats)
    r_marg = pearson(lcn, d_log)
    r_wb   = pearson(lcn_wb, d_wb)
    # Partial correlation with log2(CN) controlling for spliced BZLF1 RPM (within pop+batch)
    g_res_bz = resid(d_wb, bz_wb)
    r_part_bz = pearson(lcn_res_bz, g_res_bz)
    results.append((g, kclass, glen, st.mean(d_rpkm), st.median(d_rpkm), r_tech, r_bio, r_marg, r_wb, r_part_bz, d_rpkm, d_wb))

print('=== AIM 1.1 & 1.3: TOP 30 EBV GENES RANKED BY WITHIN-POP+BATCH CORRELATION WITH log2(CN) (n=717 clean donors) ===')
print(f'{"Gene":<13} {"Class":<12} {"Len":>5} {"MeanRPKM":>8} {"MedRPKM":>8} {"r_tech":>7} {"r_bio":>7} {"r_CN":>7} {"r_CN_wb":>8} {"r_part|BZ":>9}')
for row in sorted(results, key=lambda x: x[8], reverse=True)[:30]:
    print(f'{row[0]:<13} {row[1]:<12} {row[2]:5d} {row[3]:8.2f} {row[4]:8.2f} {row[5]:+7.3f} {row[6]:+7.3f} {row[7]:+7.4f} {row[8]:+8.4f} {row[9]:+9.4f}')

print('\n=== ALL LATENT, EBER, BART, AND B95-8 DELETION GENES (n=717 clean donors) ===')
print(f'{"Gene":<13} {"Class":<12} {"Len":>5} {"MeanRPKM":>8} {"MedRPKM":>8} {"r_tech":>7} {"r_bio":>7} {"r_CN":>7} {"r_CN_wb":>8} {"r_part|BZ":>9}')
for row in [r for r in results if r[1] in ('LATENT_MRNA', 'LATENT_BART', 'POL3_EBER', 'B958_DEL')]:
    print(f'{row[0]:<13} {row[1]:<12} {row[2]:5d} {row[3]:8.2f} {row[4]:8.2f} {row[5]:+7.3f} {row[6]:+7.3f} {row[7]:+7.4f} {row[8]:+8.4f} {row[9]:+9.4f}')

# Multiple R^2 of log2(CN) explained by Spliced BZLF1 alone vs adding the 6 Core Lytic DNA Replication Genes + BHRF1/BALF1
g_map = {r[0]: r[11] for r in results} # within-pop+batch demeaned log2 RPKM
r_bz_only = pearson(lcn_wb, bz_wb)
print(f'\n=== AIM 1.3: ABORTIVE (BZLF1) vs PRODUCTIVE LYTIC REPLICATION FORK GENES IN PREDICTING log2(CN) ===')
print(f'  Spliced BZLF1 alone within-pop+batch R^2 = {r_bz_only**2:.4f} (r = {r_bz_only:+.4f})')
for g in ['BRLF1', 'BALF5', 'BMRF1', 'BALF2', 'BBLF4', 'BSLF1', 'BBLF2_BBLF3', 'BHRF1', 'BALF1', 'BGLF4', 'BGLF5', 'BcLF1', 'BLLF1', 'BALF4']:
    r_p = pearson(lcn_res_bz, resid(g_map[g], bz_wb))
    print(f'  + {g:<12} partial r | BZLF1 = {r_p:+.4f} (t = {r_p*math.sqrt(714/(1-r_p**2)):+6.2f})')

# Aim 1.2: 3'-Coterminal Poly(A) Step-Up Flux Analysis
raw_map = {r[0]: r[10] for r in results}
clusters = [
    ('IE Zp/Rp (- strand)', ['BRLF1', 'BZLF1']),
    ('BFRF (+ strand)',     ['BFRF1A', 'BFRF1', 'BFRF2', 'BFRF3']),
    ('BMRF (+ strand)',     ['BORF2', 'BaRF1', 'BMRF1', 'BMRF2']),
    ('BGLF (- strand)',     ['BGLF5', 'BGLF4', 'BGLF3.5', 'BGLF3', 'BDLF4', 'BDLF3.5']),
]
print('\n=== AIM 1.2: 3\'-COTERMINAL POLY(A) CLUSTER READ-THROUGH vs STEP-UP PROMOTER FLUX ===')
for cname, cgenes in clusters:
    print(f'  Cluster {cname}:')
    prev = [0.0]*len(clean_donors)
    for idx, cg in enumerate(cgenes):
        cur = raw_map[cg]
        step = [max(0.0, c - p) for c, p in zip(cur, prev)]
        step_log_wb = demean(demean([math.log2(v + 0.1) for v in step], pops), bats)
        print(f'    [{idx+1}] {cg:<10}: mean raw RPKM = {st.mean(cur):7.2f} (r_CN_wb = {pearson(lcn_wb, g_map[cg]):+.4f}) | step-up flux RPKM = {st.mean(step):7.2f} (r_CN_wb = {pearson(lcn_wb, step_log_wb):+.4f})')
        prev = cur

with open(f'{BASE}/metadata/viral_gene_dissection_84.tsv', 'w') as out:
    out.write('gene\tkinetic_class\tlength_bp\tmean_rpkm\tmedian_rpkm\tr_tech\tr_bio\tr_cn_marginal\tr_cn_within_pop_batch\tr_part_given_bzlf1\n')
    for row in sorted(results, key=lambda x: x[8], reverse=True):
        out.write(f'{row[0]}\t{row[1]}\t{row[2]}\t{row[3]:.4f}\t{row[4]:.4f}\t{row[5]:.4f}\t{row[6]:.4f}\t{row[7]:.4f}\t{row[8]:.4f}\t{row[9]:.4f}\n')
