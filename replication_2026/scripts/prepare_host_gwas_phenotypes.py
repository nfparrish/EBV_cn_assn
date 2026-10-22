#!/usr/bin/env python3
"""
1. Consolidate 3,200-sample viral covariates (IR1 repeat copy number + 11 viral SNPs).
2. Run Layer A: Viral-conditioned trio heritability decomposition across 570 clean trios.
3. Build standardized residualized phenotype table for PLINK2 genome-wide & targeted host scans:
   - WGS_CN_PopAdj (n=2,449)
   - WGS_CN_ViralCond (n=2,449, conditioned on IR1 CN + 12,129:T>A + 10 viral SNPs + GC-bias)
   - WGS_IR1_PopAdj (n=2,449, specificity control)
   - RNA_BGLF5_ViralCond (n=717)
   - RNA_BZLF1_ViralCond (n=717)
   - RNA_IELatRatio_ViralCond (n=717)
   - RNA_ProdFork_ViralCond (n=717, BALF5+BcLF1 conditioned on BZLF1 + viral genetics)
   - RNA_LMP1_ViralCond (n=717)
   - RNA_Latent_ViralCond (n=717)
   - RNA_BGLF5_PopBatchAdj (n=717, unconditioned on viral genetics)
   - RNA_LMP1_PopBatchAdj (n=717, unconditioned on viral genetics)
"""
import csv
import glob
import math
import os
from collections import defaultdict
import numpy as np

try:
    from scipy import stats
except Exception:
    class _StatsFallback:
        @staticmethod
        def pearsonr(x, y):
            xa = np.asarray(x, dtype=float)
            ya = np.asarray(y, dtype=float)
            xa = xa - xa.mean()
            ya = ya - ya.mean()
            d = float(np.sqrt(np.dot(xa, xa) * np.dot(ya, ya)))
            r = float(np.dot(xa, ya) / d) if d > 0 else 0.0
            return r, 0.0
    stats = _StatsFallback()

BASE = "/work/nfp8/EBV"
META = f"{BASE}/metadata"
GWAS_DIR = f"{BASE}/host_gwas"

VIRAL_COLS = [
    "v12129_TA", "v163858_GT", "v91192_GA", "v132965_GGC",
    "v162476_TC", "v169200_CT", "v38566_GA", "v65013_CT",
    "v81920_CT", "v103742_AG", "v147913_GA"
]


def ols_residualize_z(y, X_Continuous, cat_groups_list):
    """Residualize y on categorical group lists (one-hot) and continuous covariates, then Z-standardize."""
    n = len(y)
    y_arr = np.asarray(y, dtype=float)
    cols = [np.ones(n, dtype=float)]
    for grp in cat_groups_list:
        uniq = sorted(set(grp))
        for u in uniq[1:]:
            cols.append(np.array([1.0 if g == u else 0.0 for g in grp], dtype=float))
    if X_Continuous is not None and len(X_Continuous) > 0:
        for col in X_Continuous:
            c = np.asarray(col, dtype=float)
            # mean impute any NaNs
            mask = np.isnan(c)
            if np.any(mask):
                c = c.copy()
                c[mask] = np.nanmean(c) if not np.all(mask) else 0.0
            cols.append(c)
    X = np.column_stack(cols)
    beta, _, _, _ = np.linalg.lstsq(X, y_arr, rcond=None)
    resid = y_arr - X @ beta
    sd = np.std(resid, ddof=1)
    return (resid - np.mean(resid)) / (sd if sd > 0 else 1.0), resid


def main():
    os.makedirs(GWAS_DIR, exist_ok=True)
    # 1. Consolidate viral_covars_chunks
    chunk_files = sorted(glob.glob(f"{GWAS_DIR}/viral_covars_chunks/chunk_*.tsv"))
    vcov = {}
    header = None
    for cf in chunk_files:
        with open(cf) as f:
            for idx, line in enumerate(f):
                p = line.strip().split("\t")
                if idx == 0:
                    header = p
                elif len(p) >= 5:
                    vcov[p[0]] = dict(zip(header, p))

    out_vcov = f"{META}/all_viral_covars_3200.tsv"
    if header:
        with open(out_vcov, "w") as out:
            out.write("\t".join(header) + "\n")
            for s in sorted(vcov):
                out.write("\t".join(vcov[s][h] for h in header) + "\n")
    print(f"Consolidated viral covariates for {len(vcov)} samples -> {out_vcov}")

    # Load pedigree, CN, WT fraction, GC-bias
    ped = {}
    with open(f"{META}/ped3202.txt") as f:
        for line in f:
            if line.startswith("FamilyID"):
                continue
            p = line.split()
            if len(p) >= 7:
                ped[p[1]] = {"fam": p[0], "dad": p[2], "mom": p[3], "sex": p[4], "pop": p[5], "cont": p[6]}

    cn30, wt30 = {}, {}
    with open(f"{META}/all_cn_3200.tsv") as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 10:
                cn30[p[0]] = float(p[3])
                wt30[p[0]] = float(p[9])

    gcbias = {}
    with open(f"{META}/all_gcbias_2502.tsv") as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 5:
                gcbias[p[0]] = float(p[4])

    # Load host genotype PCs from 1kgp_top_pcs.tsv and mage_eQTL_covariates.tsv
    host_pcs = {}
    if os.path.exists(f"{GWAS_DIR}/1kgp_top_pcs.tsv"):
        with open(f"{GWAS_DIR}/1kgp_top_pcs.tsv") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                s = row["sampleID"]
                host_pcs[s] = [float(row[f"PC{i}"]) for i in range(1, 11)]
    mage_pcs = {}
    if os.path.exists(f"{GWAS_DIR}/mage_eQTL_covariates.tsv"):
        with open(f"{GWAS_DIR}/mage_eQTL_covariates.tsv") as f:
            lines = [l.strip().split("\t") for l in f if l.strip()]
            samples_hdr = lines[0][1:] if len(lines[0]) == len(lines[1]) else lines[0]
            cov_dict = {l[0]: [float(x) for x in l[1:]] for l in lines[1:] if l[0].startswith("PC")}
            for idx, s in enumerate(samples_hdr):
                if "PC1" in cov_dict:
                    mage_pcs[s] = [cov_dict[f"PC{i}"][idx] for i in range(1, 16) if f"PC{i}" in cov_dict]

    # Helper to get viral feature vector per sample
    def get_ir1(s):
        v = vcov.get(s, {}).get("ir1_copies", "NA")
        return float(v) if v != "NA" and float(v) > 0 else np.nan

    def get_vsnp(s, col):
        v = vcov.get(s, {}).get(col, "NA")
        return float(v) if v != "NA" else np.nan

    # =========================================================================
    # LAYER A: Viral-Conditioned Trio Heritability (n=570 clean trios)
    # =========================================================================
    clean_3200 = sorted([s for s in cn30 if wt30.get(s, 0) <= 0.01 and cn30[s] > 0 and s in ped and s in vcov])
    y_3200 = np.array([math.log2(cn30[s]) for s in clean_3200])
    pops_3200 = [ped[s]["pop"] for s in clean_3200]
    sex_3200 = [ped[s]["sex"] for s in clean_3200]
    ir1_3200 = [get_ir1(s) for s in clean_3200]
    v12129_3200 = [get_vsnp(s, "v12129_TA") for s in clean_3200]
    v163858_3200 = [get_vsnp(s, "v163858_GT") for s in clean_3200]
    all_vsnps_3200 = [[get_vsnp(s, col) for s in clean_3200] for col in VIRAL_COLS]

    # Build sample -> residualized trait maps for 5 nested models
    _, res_m1 = ols_residualize_z(y_3200, [], [pops_3200])
    _, res_m2 = ols_residualize_z(y_3200, [ir1_3200], [pops_3200])
    _, res_m3 = ols_residualize_z(y_3200, [ir1_3200, v12129_3200, v163858_3200], [pops_3200])
    _, res_m4 = ols_residualize_z(y_3200, [ir1_3200] + all_vsnps_3200, [pops_3200, sex_3200])
    _, res_ir1_w = ols_residualize_z(np.nan_to_num(ir1_3200, nan=np.nanmean(ir1_3200)), [], [pops_3200])

    idx_map = {s: i for i, s in enumerate(clean_3200)}
    trios = []
    for s, info in ped.items():
        fa, mo, pop = info["dad"], info["mom"], info["pop"]
        if fa != "0" and mo != "0" and s in idx_map and fa in idx_map and mo in idx_map:
            trios.append((s, fa, mo, pop))

    def eval_trio_model(name, val_arr):
        fa_v = np.array([val_arr[idx_map[fa]] for _, fa, _, _ in trios])
        mo_v = np.array([val_arr[idx_map[mo]] for _, _, mo, _ in trios])
        ch_v = np.array([val_arr[idx_map[s]] for s, _, _, _ in trios])
        mid_v = 0.5 * (fa_v + mo_v)
        r_sp = stats.pearsonr(np.concatenate([fa_v, mo_v]), np.concatenate([ch_v, ch_v]))[0]
        r_pp = stats.pearsonr(fa_v, mo_v)[0]
        slope_mp = np.cov(ch_v, mid_v, ddof=1)[0, 1] / np.var(mid_v, ddof=1)
        h2_adj = 2.0 * (r_sp - r_pp)
        # Bootstrap SE (1,000 resamples)
        rng = np.random.default_rng(42)
        n_t = len(trios)
        boot_h2 = []
        for _ in range(1000):
            b = rng.integers(0, n_t, size=n_t)
            rsp_b = stats.pearsonr(np.concatenate([fa_v[b], mo_v[b]]), np.concatenate([ch_v[b], ch_v[b]]))[0]
            rpp_b = stats.pearsonr(fa_v[b], mo_v[b])[0]
            boot_h2.append(2.0 * (rsp_b - rpp_b))
        se_h2 = float(np.std(boot_h2, ddof=1))
        ci_lo, ci_hi = float(np.percentile(boot_h2, 2.5)), float(np.percentile(boot_h2, 97.5))
        return {
            "model": name, "n_trios": n_t, "r_parent_child": r_sp, "r_spouse": r_pp,
            "midparent_slope": slope_mp, "h2_spouse_corrected": h2_adj,
            "h2_se": se_h2, "h2_ci95_lo": ci_lo, "h2_ci95_hi": ci_hi,
            "resid_sd": float(np.std(val_arr, ddof=1))
        }

    trio_results = [
        eval_trio_model("M0_Raw_Pooled_log2CN", y_3200),
        eval_trio_model("M1_Within_Pop_log2CN", res_m1),
        eval_trio_model("M2_Within_Pop_plus_IR1_CN", res_m2),
        eval_trio_model("M3_Within_Pop_plus_IR1_and_12129_163858", res_m3),
        eval_trio_model("M4_Within_Pop_plus_Full_Viral_Genetics_11SNPs_IR1", res_m4),
        eval_trio_model("Viral_IR1_Repeat_Count_Within_Pop_Control", res_ir1_w),
    ]

    trio_out = f"{META}/trio_heritability_viral_conditioned.tsv"
    with open(trio_out, "w") as out:
        cols = list(trio_results[0].keys())
        out.write("\t".join(cols) + "\n")
        for r in trio_results:
            out.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
            print(f"  [Trio {r['model']:<48s}] r_sp={r['r_parent_child']:+.4f}, r_pp={r['r_spouse']:+.4f}, h2={r['h2_spouse_corrected']:+.4f} (95% CI [{r['h2_ci95_lo']:+.3f}, {r['h2_ci95_hi']:+.3f}]), SD={r['resid_sd']:.3f}")

    # =========================================================================
    # BUILD RESIDUALIZED PHENOTYPES FOR WGS (n=2,449) & MAGE RNA (n=717)
    # =========================================================================
    unrel = set()
    with open(f"{META}/ebv_cn_1kGP_2504.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            unrel.add(row["sample"])
    clean_2449 = sorted([s for s in unrel if s in cn30 and s in ped and wt30.get(s, 0) <= 0.01])

    y_wgs = np.array([math.log2(cn30[s]) for s in clean_2449])
    pops_wgs = [ped[s]["pop"] for s in clean_2449]
    sex_wgs = [ped[s]["sex"] for s in clean_2449]
    pcs_wgs = [[host_pcs.get(s, [0.0]*10)[k] for s in clean_2449] for k in range(10)]
    ir1_wgs = [get_ir1(s) for s in clean_2449]
    gc_wgs = [gcbias.get(s, np.nan) for s in clean_2449]
    vsnps_wgs = [[get_vsnp(s, col) for s in clean_2449] for col in VIRAL_COLS]

    z_wgs_pop, _ = ols_residualize_z(y_wgs, pcs_wgs, [pops_wgs, sex_wgs])
    z_wgs_vcond, _ = ols_residualize_z(y_wgs, pcs_wgs + [ir1_wgs, gc_wgs] + vsnps_wgs, [pops_wgs, sex_wgs])
    z_wgs_ir1, _ = ols_residualize_z(np.nan_to_num(ir1_wgs, nan=np.nanmean(ir1_wgs)), pcs_wgs, [pops_wgs, sex_wgs])

    # MAGE 717 clean donors
    by_donor_rna = defaultdict(list)
    with open(f"{META}/all_mage_ebv_rna_779.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by_donor_rna[row["donor"]].append(row)
    by_donor_genes = defaultdict(lambda: defaultdict(list))
    with open(f"{META}/all_mage_ebv_genes_779.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            by_donor_genes[row["donor"]][row["gene"]].append(float(row["rpkm"]))

    clean_717 = sorted([d for d in by_donor_rna if d in cn30 and d in ped and wt30.get(d, 0) <= 0.01])
    pops_m = [ped[d]["pop"] for d in clean_717]
    sex_m = [ped[d]["sex"] for d in clean_717]
    batch_m = [str(by_donor_rna[d][0]["batch"]) for d in clean_717]
    def get_pcs10(s):
        v = mage_pcs.get(s, host_pcs.get(s, []))
        return (list(v) + [0.0] * 10)[:10]

    pcs_m = [[get_pcs10(d)[k] for d in clean_717] for k in range(10)]
    ir1_m = [get_ir1(d) for d in clean_717]
    gc_m = [gcbias.get(d, np.nan) for d in clean_717]
    vsnps_m = [[get_vsnp(d, col) for d in clean_717] for col in VIRAL_COLS]

    y_bglf5 = np.array([math.log2(np.mean(by_donor_genes[d]["BGLF5"]) + 0.01) for d in clean_717])
    y_bzlf1 = np.array([math.log2(np.mean([float(r["bzlf1_spliced_reads"]) / float(r["total_reads"]) * 1e6 for r in by_donor_rna[d]]) + 0.001) for d in clean_717])
    y_ielat = np.array([np.mean([float(r["lytic_ie_over_latent_log2"]) for r in by_donor_rna[d]]) for d in clean_717])
    y_balf5 = np.array([math.log2(np.mean(by_donor_genes[d]["BALF5"]) + 0.01) for d in clean_717])
    y_bclf1 = np.array([math.log2(np.mean(by_donor_genes[d]["BcLF1"]) + 0.01) for d in clean_717])
    y_prodfork = 0.5 * (y_balf5 + y_bclf1)
    y_lmp1 = np.array([math.log2(np.mean(by_donor_genes[d]["LMP-1"]) + 0.01) for d in clean_717])
    y_latent = np.array([math.log2(np.mean([float(r["latent_mrna_rpkm"]) for r in by_donor_rna[d]]) + 0.01) for d in clean_717])

    vcov_list_m = pcs_m + [ir1_m, gc_m] + vsnps_m
    z_bglf5_vc, _ = ols_residualize_z(y_bglf5, vcov_list_m, [pops_m, batch_m, sex_m])
    z_bzlf1_vc, _ = ols_residualize_z(y_bzlf1, vcov_list_m, [pops_m, batch_m, sex_m])
    z_ielat_vc, _ = ols_residualize_z(y_ielat, vcov_list_m, [pops_m, batch_m, sex_m])
    z_pfork_vc, _ = ols_residualize_z(y_prodfork, vcov_list_m + [y_bzlf1], [pops_m, batch_m, sex_m])
    z_lmp1_vc, _ = ols_residualize_z(y_lmp1, vcov_list_m, [pops_m, batch_m, sex_m])
    z_lat_vc, _ = ols_residualize_z(y_latent, vcov_list_m, [pops_m, batch_m, sex_m])

    z_bglf5_uncond, _ = ols_residualize_z(y_bglf5, pcs_m, [pops_m, batch_m, sex_m])
    z_lmp1_uncond, _ = ols_residualize_z(y_lmp1, pcs_m, [pops_m, batch_m, sex_m])

    # Write combined PLINK2 phenotype table (#IID + 11 phenotypes)
    wgs_idx = {s: i for i, s in enumerate(clean_2449)}
    rna_idx = {d: i for i, d in enumerate(clean_717)}
    all_samples = sorted(set(clean_2449) | set(clean_717))

    pheno_names = [
        "WGS_CN_PopAdj", "WGS_CN_ViralCond", "WGS_IR1_PopAdj",
        "RNA_BGLF5_ViralCond", "RNA_BZLF1_ViralCond", "RNA_IELatRatio_ViralCond",
        "RNA_ProdFork_ViralCond", "RNA_LMP1_ViralCond", "RNA_Latent_ViralCond",
        "RNA_BGLF5_Uncond", "RNA_LMP1_Uncond"
    ]
    plink_pheno_path = f"{GWAS_DIR}/ebv_host_phenotypes.plink.tsv"
    with open(plink_pheno_path, "w") as out:
        out.write("#IID\t" + "\t".join(pheno_names) + "\n")
        for s in all_samples:
            vals = []
            if s in wgs_idx:
                i = wgs_idx[s]
                vals.extend([f"{z_wgs_pop[i]:.5f}", f"{z_wgs_vcond[i]:.5f}", f"{z_wgs_ir1[i]:.5f}"])
            else:
                vals.extend(["NA", "NA", "NA"])
            if s in rna_idx:
                j = rna_idx[s]
                vals.extend([
                    f"{z_bglf5_vc[j]:.5f}", f"{z_bzlf1_vc[j]:.5f}", f"{z_ielat_vc[j]:.5f}",
                    f"{z_pfork_vc[j]:.5f}", f"{z_lmp1_vc[j]:.5f}", f"{z_lat_vc[j]:.5f}",
                    f"{z_bglf5_uncond[j]:.5f}", f"{z_lmp1_uncond[j]:.5f}"
                ])
            else:
                vals.extend(["NA"] * 8)
            out.write(s + "\t" + "\t".join(vals) + "\n")
    print(f"Wrote PLINK2 phenotype matrix ({len(all_samples)} samples x {len(pheno_names)} traits) -> {plink_pheno_path}")


if __name__ == "__main__":
    main()
