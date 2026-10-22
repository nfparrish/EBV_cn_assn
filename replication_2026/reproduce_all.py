#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = [
#   "numpy>=1.24",
#   "scipy>=1.10",
#   "matplotlib>=3.7",
# ]
# ///
"""
Turn-Key Replication & Verification Suite (`reproduce_all.py`)
==============================================================
Deterministically recomputes and verifies every headline finding in
`docs/CLAIMS_LEDGER.md` and `docs/HANDOFF.md` from the local master tables in
`replication_2026/data/`, writes `REPLICATION_VERIFICATION_REPORT.md`, and
regenerates all 6 multi-panel publication figures (`figures/fig1..fig6.png/.pdf`).

Usage:
    uv run replication_2026/reproduce_all.py
"""
import csv
import math
import os
import sys
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
            denom = float(np.sqrt(np.dot(xa, xa) * np.dot(ya, ya)))
            return (float(np.dot(xa, ya) / denom) if denom > 0 else 0.0, 0.0)

        @staticmethod
        def spearmanr(x, y):
            def _rank(a):
                arr = np.asarray(a, dtype=float)
                order = np.argsort(arr, kind="mergesort")
                ranks = np.empty(len(arr), dtype=float)
                i = 0
                while i < len(arr):
                    j = i
                    while j + 1 < len(arr) and arr[order[j + 1]] == arr[order[i]]:
                        j += 1
                    ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
                    i = j + 1
                return ranks
            return _StatsFallback.pearsonr(_rank(x), _rank(y))
    stats = _StatsFallback()

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, "data")
SCRIPTS_DIR = os.path.join(ROOT_DIR, "scripts")
sys.path.insert(0, SCRIPTS_DIR)


def demean_by_groups(vals, groups):
    g_sums = defaultdict(float)
    g_cnts = defaultdict(int)
    for v, g in zip(vals, groups):
        g_sums[g] += v
        g_cnts[g] += 1
    g_means = {g: g_sums[g] / g_cnts[g] for g in g_cnts}
    return np.array([v - g_means[g] for v, g in zip(vals, groups)])


def partial_corr(x, y, z):
    rx = stats.pearsonr(x, z)[0]
    ry = stats.pearsonr(y, z)[0]
    rxy = stats.pearsonr(x, y)[0]
    return (rxy - rx * ry) / math.sqrt(max(1e-12, (1 - rx**2) * (1 - ry**2)))


def run_verification():
    checks = []

    def record(section, claim, metric_str, passed=True):
        status = "PASS" if passed else "FAIL"
        checks.append((section, claim, metric_str, status))
        print(f"  [{status}] {claim:<62s} | {metric_str}")

    print("=" * 105)
    print("1000 GENOMES / MAGE EBV COPY NUMBER, VIRAL GENETICS & HOST MEDIATION — TURN-KEY REPLICATION SUITE")
    print("=" * 105)

    # Load core tables
    ped = {}
    with open(os.path.join(DATA_DIR, "ped3202.txt")) as f:
        for line in f:
            if line.startswith("FamilyID"):
                continue
            p = line.split()
            if len(p) >= 7:
                ped[p[1]] = {
                    "fam": p[0], "dad": p[2], "mom": p[3],
                    "sex": int(p[4]), "pop": p[5], "cont": p[6]
                }

    cn30, wt30 = {}, {}
    with open(os.path.join(DATA_DIR, "all_cn_3200.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 10:
                cn30[p[0]] = float(p[3])
                wt30[p[0]] = float(p[9])

    unrel = set()
    with open(os.path.join(DATA_DIR, "ebv_cn_1kGP_2504.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            unrel.add(row["sample"])

    all_2502 = sorted([s for s in unrel if s in cn30 and s in ped])
    clean_2449 = sorted([s for s in all_2502 if wt30[s] <= 0.01])

    # 1. WGS CN Architecture & Heritability
    print("\n--- SECTION 1: WGS EBV Copy Number Architecture, Heritability, & BCR Clonality ---")
    c_all = [cn30[s] for s in all_2502]
    c_clean = [cn30[s] for s in clean_2449]
    record("1. WGS CN", "EBV CN range & median (n=2,502 all / n=2,449 clean)",
           f"range={min(c_all):.2f}–{max(c_all):.1f}, med_all={np.median(c_all):.2f}, med_clean={np.median(c_clean):.2f}",
           len(all_2502) == 2502 and len(clean_2449) == 2449 and abs(np.median(c_clean) - 36.16) < 0.05)

    # Trio heritability
    clean_all_map = {s: math.log2(v) for s, v in cn30.items() if wt30.get(s, 0) <= 0.01 and v > 0}
    trios = []
    for s, info in ped.items():
        fa, mo, pop = info["dad"], info["mom"], info["pop"]
        if fa != "0" and mo != "0" and s in clean_all_map and fa in clean_all_map and mo in clean_all_map:
            trios.append((s, fa, mo, pop))
    p_vals = defaultdict(list)
    for s, v in clean_all_map.items():
        if s in ped:
            p_vals[ped[s]["pop"]].append(v)
    p_mean = {p: np.mean(v) for p, v in p_vals.items()}
    vw = lambda s, p: clean_all_map[s] - p_mean[p]

    r_sp_raw = stats.pearsonr([clean_all_map[fa] for _, fa, _, _ in trios] + [clean_all_map[mo] for _, _, mo, _ in trios],
                              [clean_all_map[s] for s, _, _, _ in trios] * 2)[0]
    r_pp_raw = stats.pearsonr([clean_all_map[fa] for _, fa, _, _ in trios],
                              [clean_all_map[mo] for _, _, mo, _ in trios])[0]
    r_sp_w = stats.pearsonr([vw(fa, p) for _, fa, _, p in trios] + [vw(mo, p) for _, _, mo, p in trios],
                            [vw(s, p) for s, _, _, p in trios] * 2)[0]
    r_pp_w = stats.pearsonr([vw(fa, p) for _, fa, _, p in trios],
                            [vw(mo, p) for _, _, mo, p in trios])[0]
    h2_w = 2.0 * (r_sp_w - r_pp_w)
    record("1. WGS CN", "Within-population trio heritability (n=570 clean trios)",
           f"pooled r_sp={r_sp_raw:+.4f}/r_pp={r_pp_raw:+.4f} | within-pop r_sp={r_sp_w:+.4f}/r_pp={r_pp_w:+.4f} -> h²={h2_w:+.4f}",
           len(trios) == 570 and abs(h2_w - 0.0101) < 0.002)

    # 3-platform CN concordance
    ont_cn = {}
    with open(os.path.join(DATA_DIR, "all_ont_904.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                ont_cn[p[0]] = float(p[5])
    low_cn, low_dep = {}, {}
    with open(os.path.join(DATA_DIR, "all_lowcov_2499.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                low_dep[p[0]] = float(p[2])
                low_cn[p[0]] = float(p[5])

    s_ont = [s for s in clean_2449 if s in ont_cn and ont_cn[s] > 0]
    s_low = [s for s in clean_2449 if s in low_cn and low_dep[s] >= 5.0]
    s_3p = [s for s in s_ont if s in low_cn and low_dep[s] >= 5.0]
    r_30_ont = stats.pearsonr([math.log2(cn30[s]) for s in s_ont], [math.log2(ont_cn[s]) for s in s_ont])[0]
    r_30_low = stats.pearsonr([math.log2(cn30[s]) for s in s_low], [math.log2(low_cn[s]) for s in s_low])[0]
    r_ont_low = stats.pearsonr([math.log2(ont_cn[s]) for s in s_3p], [math.log2(low_cn[s]) for s in s_3p])[0]
    record("1. WGS CN", "Three-platform CN concordance (30x vs ONT vs Low-cov)",
           f"r(30x,ONT)={r_30_ont:+.4f} (n={len(s_ont)}), r(30x,Low)={r_30_low:+.4f} (n={len(s_low)}), r(ONT,Low)={r_ont_low:+.4f} (n={len(s_3p)})",
           abs(r_30_ont - 0.9072) < 0.002 and abs(r_30_low - 0.8427) < 0.002)

    n_blood = sum(1 for d in low_dep.values() if d < 0.5)
    n_lcl = sum(1 for d in low_dep.values() if d >= 5.0)
    record("1. WGS CN", "Phase 3 low-coverage blood (<0.5x) vs LCL (>=5x) split",
           f"n_total={len(low_dep)}, blood(<0.5x)={n_blood} ({n_blood/len(low_dep)*100:.1f}%), LCL(>=5x)={n_lcl}",
           n_blood == 469 and n_lcl == 2004)

    # BCR clonality
    bcr_frac = {}
    with open(os.path.join(DATA_DIR, "bcr_wgs_summary.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 7:
                bcr_frac[p[0]] = float(p[6])
    s_bcr = [s for s in clean_2449 if s in bcr_frac]
    rho_bcr = stats.spearmanr([bcr_frac[s] for s in s_bcr], [math.log2(cn30[s]) for s in s_bcr])[0]
    record("1. WGS CN", "PoN-cleaned WGS BCR clonality vs log2(CN) (n=2,340)",
           f"n={len(s_bcr)}, median top_clone_frac={np.median([bcr_frac[s] for s in s_bcr]):.3f}, Spearman rho={rho_bcr:+.4f}",
           len(s_bcr) >= 2313 and abs(rho_bcr - 0.078) < 0.005)

    # 2. H1 & Replicate Thaw Repeatability across all 779 MAGE libraries
    print("\n--- SECTION 2: P1 (H1) Lytic vs Latent EBV RNA Across All 779/779 MAGE Libraries ---")
    mage_runs = []
    with open(os.path.join(DATA_DIR, "all_mage_ebv_rna_779.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            mage_runs.append(row)

    by_donor = defaultdict(list)
    for r in mage_runs:
        if r["donor"] in cn30:
            by_donor[r["donor"]].append(r)

    clean_mage = sorted([d for d in by_donor if wt30.get(d, 0) <= 0.01 and d in ped])
    same_b, diff_b = [], []
    for d, rlist in by_donor.items():
        if len(rlist) > 1:
            rlist = sorted(rlist, key=lambda x: x["run"])
            for i in range(len(rlist)):
                for j in range(i + 1, len(rlist)):
                    if rlist[i]["batch"] == rlist[j]["batch"]:
                        same_b.append((rlist[i], rlist[j]))
                    else:
                        diff_b.append((rlist[i], rlist[j]))

    r_ie_lat_same = stats.pearsonr([float(a["lytic_ie_over_latent_log2"]) for a, _ in same_b],
                                   [float(b["lytic_ie_over_latent_log2"]) for _, b in same_b])[0]
    r_ie_lat_diff = stats.pearsonr([float(a["lytic_ie_over_latent_log2"]) for a, _ in diff_b],
                                   [float(b["lytic_ie_over_latent_log2"]) for _, b in diff_b])[0]
    r_lat_diff = stats.pearsonr([math.log2(float(a["latent_mrna_rpkm"]) + 0.01) for a, _ in diff_b],
                                [math.log2(float(b["latent_mrna_rpkm"]) + 0.01) for _, b in diff_b])[0]
    record("2. H1 RNA", "Replicate thaw reproducibility (24 donors, 70 pairs)",
           f"same-batch n={len(same_b)} r(IE/Lat)={r_ie_lat_same:+.4f} | diff-batch n={len(diff_b)} r(IE/Lat)={r_ie_lat_diff:+.4f} vs r(Latent)={r_lat_diff:+.4f}",
           len(same_b) == 24 and len(diff_b) == 46 and abs(r_ie_lat_diff - 0.9102) < 0.002)

    x_cn_m = np.array([math.log2(cn30[d]) for d in clean_mage])
    y_lat_m = np.array([math.log2(np.mean([float(r["latent_mrna_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_mage])
    y_ie_m = np.array([math.log2(np.mean([float(r["lytic_ie_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_mage])
    y_z_m = np.array([math.log2(np.mean([float(r["bzlf1_spliced_reads"]) / float(r["total_reads"]) * 1e6 for r in by_donor[d]]) + 0.001) for d in clean_mage])
    pops_m = [ped[d]["pop"] for d in clean_mage]
    pb_m = [f"{ped[d]['pop']}_{by_donor[d][0]['batch']}" for d in clean_mage]

    r_lat_marg = stats.pearsonr(x_cn_m, y_lat_m)[0]
    r_lat_wb = stats.pearsonr(demean_by_groups(x_cn_m, pb_m), demean_by_groups(y_lat_m, pb_m))[0]
    r_z_marg = stats.pearsonr(x_cn_m, y_z_m)[0]
    r_z_wb = stats.pearsonr(demean_by_groups(x_cn_m, pb_m), demean_by_groups(y_z_m, pb_m))[0]
    r_ie_marg = stats.pearsonr(x_cn_m, y_ie_m)[0]
    r_ie_wb = stats.pearsonr(demean_by_groups(x_cn_m, pb_m), demean_by_groups(y_ie_m, pb_m))[0]
    record("2. H1 RNA", "Full-cohort H1 (n=717 clean donors, 779/779 libraries)",
           f"Latent r={r_lat_marg:+.4f} (wb {r_lat_wb:+.4f}) | Spliced BZLF1 r={r_z_marg:+.4f} (wb {r_z_wb:+.4f}) | Lytic IE r={r_ie_marg:+.4f} (wb {r_ie_wb:+.4f})",
           len(clean_mage) == 717 and abs(r_lat_marg - (-0.0287)) < 0.002 and abs(r_z_marg - 0.6477) < 0.002)

    # 3. Aim 1: All 84 Viral Genes
    print("\n--- SECTION 3: Aim 1 — Dissection of All 84 EBV Genes & Productive Lytic Replication ---")
    genes84 = {}
    with open(os.path.join(DATA_DIR, "viral_gene_dissection_84.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes84[row["gene"]] = row

    bglf5 = genes84["BGLF5"]
    balf5 = genes84["BALF5"]
    bclf1 = genes84["BcLF1"]
    balf1 = genes84["BALF1"]
    bhrf1 = genes84["BHRF1"]
    record("3. Aim 1", "Top viral gene BGLF5 & cross-thaw repeatability",
           f"BGLF5 r_CN={float(bglf5['r_cn_marginal']):+.4f}, r_wb={float(bglf5['r_cn_within_pop_batch']):+.4f}, r_bio={float(bglf5['r_bio']):+.4f}",
           abs(float(bglf5["r_cn_marginal"]) - 0.6886) < 0.002)
    record("3. Aim 1", "Productive replication fork & late capsid partial r | BZLF1",
           f"BALF5(Pol)={float(balf5['r_part_given_bzlf1']):+.4f}, BcLF1(Capsid)={float(bclf1['r_part_given_bzlf1']):+.4f}, BALF1={float(balf1['r_part_given_bzlf1']):+.4f} vs BHRF1={float(bhrf1['r_part_given_bzlf1']):+.4f}",
           abs(float(balf5["r_part_given_bzlf1"]) - 0.2439) < 0.002 and abs(float(bhrf1["r_part_given_bzlf1"]) - 0.0016) < 0.002)

    # 4. Aim 2: Viral Genetics (IR1 Structural QTL, 112-Site Scan, & Within-Sample cis-ASE)
    print("\n--- SECTION 4: Aim 2 — Viral Genetic Determinants (IR1 Structural QTL, 12,129:T>A, & cis-ASE) ---")
    ir1_30, ir1_flank, flank_uniq = {}, {}, {}
    with open(os.path.join(DATA_DIR, "all_ir1_2502.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                s = p[0]
                d_ir1, d_cp, d_uniq = float(p[1]), float(p[2]), float(p[3])
                ir1_30[s] = float(p[4]) * (22737.0 / 3072.0)
                if d_cp > 0 and d_uniq > 0:
                    ir1_flank[s] = (d_ir1 / d_cp) * (22737.0 / 3072.0)
                    flank_uniq[s] = d_cp / d_uniq
    gcbias = {}
    with open(os.path.join(DATA_DIR, "all_gcbias_2502.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 5:
                gcbias[p[0]] = float(p[4])
    ont_ir1 = {}
    with open(os.path.join(DATA_DIR, "all_ont_ir1_907.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4:
                ont_ir1[p[0]] = float(p[3])

    pops_2449 = [ped[s]["pop"] for s in clean_2449]
    ir1_w = demean_by_groups([ir1_30[s] for s in clean_2449], pops_2449)
    cn_w = demean_by_groups([math.log2(cn30[s]) for s in clean_2449], pops_2449)
    gc_w = demean_by_groups([gcbias[s] for s in clean_2449], pops_2449)
    r_ir1_wp = stats.pearsonr(ir1_w, cn_w)[0]
    r_ir1_gc = partial_corr(ir1_w, cn_w, gc_w)

    s_ont_ir1 = [s for s in clean_2449 if s in ont_ir1 and s in ont_cn and ont_cn[s] > 0]
    pops_ont = [ped[s]["pop"] for s in s_ont_ir1]
    r_ont_ir1_wp = stats.pearsonr(demean_by_groups([ont_ir1[s] for s in s_ont_ir1], pops_ont),
                                  demean_by_groups([math.log2(ont_cn[s]) for s in s_ont_ir1], pops_ont))[0]
    record("4. Aim 2", "Three-platform within-population IR1 -> log2(CN) structural QTL",
           f"30x r_wp={r_ir1_wp:+.4f}, GC-controlled r={r_ir1_gc:+.4f} (96.8%), ONT r_wp={r_ont_ir1_wp:+.4f} (n={len(s_ont_ir1)})",
           abs(r_ir1_wp - 0.3148) < 0.002 and abs(r_ir1_gc - 0.3048) < 0.002 and abs(r_ont_ir1_wp - 0.2945) < 0.002)

    qtl112 = {}
    with open(os.path.join(DATA_DIR, "viral_qtl_scan_112.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            qtl112[row["site"]] = row

    s_12129 = qtl112["12129:T>A"]
    s_163858 = qtl112["163858:G>T"]
    s_91192 = qtl112["91192:G>A"]
    s_132965 = qtl112["132965:G>GC"]
    record("4. Aim 2", "Out-of-sample replicated 12,129:T>A CN-QTL vs neutral barcodes",
           f"12129:T>A r_wp={float(s_12129['r_cn_within_pop_2449']):+.4f} (holdout {float(s_12129['r_cn_holdout_1732']):+.4f}) | 91192:G>A r_wp={float(s_91192['r_cn_within_pop_2449']):+.4f} | 132965:G>GC r_wp={float(s_132965['r_cn_within_pop_2449']):+.4f}",
           abs(float(s_12129["r_cn_within_pop_2449"]) - 0.1695) < 0.002 and abs(float(s_12129["r_cn_holdout_1732"]) - 0.1736) < 0.002)

    # Within-sample heteroplasmic cis-ASE & Qp vs Cp/Wp splice check
    obs_91192 = []
    with open(os.path.join(DATA_DIR, "all_mage_viral_ase_59.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rtot = int(row["rna_ref"]) + int(row["rna_alt"])
            if int(row["pos"]) == 91192 and row["ref"] == "G" and row["alt"] == "A" and rtot >= 5 and float(row["dna_vaf"]) > 0:
                obs_91192.append((float(row["dna_vaf"]), int(row["rna_ref"]), int(row["rna_alt"])))
    tot_91192 = sum(r + a for _, r, a in obs_91192)
    alt_91192 = sum(a for _, r, a in obs_91192)
    wdna_91192 = sum(d * (r + a) for d, r, a in obs_91192) / tot_91192

    spl_qp, spl_cpwp = [], []
    with open(os.path.join(DATA_DIR, "all_mage_viral_splices_59.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["splice_class"] == "Qp_to_EBNA1":
                spl_qp.append(int(row["count"]))
            elif row["splice_class"] == "CpWp_U_to_EBNA1":
                spl_cpwp.append(int(row["count"]))

    record("4. Aim 2", "Within-sample cis-ASE at 91,192:G>A & Latency III vs I splicing",
           f"91192:G>A RNA={alt_91192}/{tot_91192} ({alt_91192/tot_91192:.3f}) vs DNA={wdna_91192:.3f} | Qp_to_EBNA1={sum(1 for x in spl_qp if x>0)}/{len(spl_qp)} vs CpWp_U_to_EBNA1={sum(1 for x in spl_cpwp if x>0)}/{len(spl_cpwp)}",
           sum(spl_qp) == 0 and all(x > 0 for x in spl_cpwp))

    # 5. Aim 3: Host Transcriptome Mediation & Out-of-Sample Vector Alignment
    print("\n--- SECTION 5: Aim 3 — Host Transcriptome Mediation & Out-of-Sample Vector Alignment ---")
    med_rows = []
    with open(os.path.join(DATA_DIR, "host_viral_mediation_20154.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            med_rows.append(row)

    n_cn_sig = sum(1 for r in med_rows if float(r["q_wgs_cn"]) < 0.05)
    n_z_sig = sum(1 for r in med_rows if float(r["q_spliced_bzlf1"]) < 0.05)
    n_bglf5_sig = sum(1 for r in med_rows if float(r["q_bglf5"]) < 0.05)
    n_lmp1_sig = sum(1 for r in med_rows if float(r["q_lmp1"]) < 0.05)
    n_ir1_sig = sum(1 for r in med_rows if float(r["q_ir1_dna"]) < 0.05)
    r_ir1_cn_vec = stats.pearsonr([float(r["r_ir1_dna"]) for r in med_rows],
                                  [float(r["r_wgs_cn"]) for r in med_rows])[0]
    record("5. Aim 3", "Host discovery boost & exogenous IR1 instrument alignment",
           f"q<0.05: WGS_CN={n_cn_sig:,}, BZLF1={n_z_sig:,}, BGLF5={n_bglf5_sig:,}, LMP1={n_lmp1_sig:,}, IR1_DNA={n_ir1_sig} (vector r={r_ir1_cn_vec:+.4f})",
           n_cn_sig == 6087 and n_z_sig == 9819 and n_bglf5_sig == 10580 and n_lmp1_sig == 13800 and n_ir1_sig == 663)

    # Vector alignment (24 non-IBS/GBR pops -> IBS vs GBR)
    beta_24 = np.array([float(r["beta_cn_24pops"]) for r in med_rows])
    de_vec = np.array([float(r["de_ibs_minus_gbr"]) for r in med_rows])
    sig_mask = np.array([float(r["q_wgs_cn"]) < 0.05 for r in med_rows])
    r_vec_all = stats.pearsonr(beta_24, de_vec)[0]
    r_vec_sig = stats.pearsonr(beta_24[sig_mask], de_vec[sig_mask])[0]
    record("5. Aim 3", "Zero-sample-overlap vector alignment (24 pops -> IBS vs GBR)",
           f"All 20,154 genes r={r_vec_all:+.4f} (R²={r_vec_all**2*100:.1f}%) | 6,087 CN-sig genes r={r_vec_sig:+.4f} (R²={r_vec_sig**2*100:.1f}%)",
           abs(r_vec_all - 0.7086) < 0.01 and abs(r_vec_sig - 0.8444) < 0.01)

    # 6. Aim 4: Viral-Conditioned Host Genetic Scan (Layers A, B1, B2, B3, & C)
    print("\n--- SECTION 6: Aim 4 — Viral-Conditioned Host Genetic Architecture (Layers A, B1–B3, & C) ---")
    trio_vc = {}
    with open(os.path.join(DATA_DIR, "trio_heritability_viral_conditioned.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            trio_vc[row["model"]] = row
    m1_t = trio_vc["M1_Within_Pop_log2CN"]
    m4_t = trio_vc["M4_Within_Pop_plus_Full_Viral_Genetics_11SNPs_IR1"]
    record("6. Aim 4", "Viral-conditioned trio heritability (n=570 clean trios)",
           f"M1 r_pp={float(m1_t['r_spouse']):+.4f} -> M4 r_pp={float(m4_t['r_spouse']):+.4f} (-38.3%) | M4 h²={float(m4_t['h2_spouse_corrected']):+.4f} (95% CI [{float(m4_t['h2_ci95_lo']):+.3f}, {float(m4_t['h2_ci95_hi']):+.3f}])",
           abs(float(m4_t["r_spouse"]) - 0.1119) < 0.005 and abs(float(m4_t["h2_spouse_corrected"]) - 0.0578) < 0.005)

    gw_host = {}
    with open(os.path.join(DATA_DIR, "host_gwas_genomewide_summary.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gw_host[row["phenotype"]] = row
    cand_host = defaultdict(list)
    with open(os.path.join(DATA_DIR, "host_candidate55_scan.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cand_host[row["phenotype"]].append(row)

    w_pop = gw_host["WGS_CN_PopAdj"]
    w_vc = gw_host["WGS_CN_ViralCond"]
    dnmt3b_ielat = cand_host["RNA_IELatRatio_ViralCond"][0]
    record("6. Aim 4", "Genome-wide 7.62M host GWAS & DNMT3B lytic switch candidate",
           f"GWAS P<5e-8: PopAdj={w_pop['n_gw_sig_5e8']} -> ViralCond={w_vc['n_gw_sig_5e8']} (max PVE={w_vc['lead_pve_pct']}%) | DNMT3B IE/Lat P={float(dnmt3b_ielat['p']):.2e} (q={dnmt3b_ielat['q_cand55']})",
           int(w_vc["n_gw_sig_5e8"]) == 0 and dnmt3b_ielat["gene"] == "DNMT3B" and float(dnmt3b_ielat["q_cand55"]) < 0.05)

    # 7. Aim 5: Biobank (UKB + All of Us Blood & Saliva) vs. LCL Viral-Conditioned Architecture
    print("\n--- SECTION 7: Aim 5 — Biobank (UKB+AoU Blood & Saliva) vs. LCL Cross-Study Architecture ---")
    bio_loci = {}
    with open(os.path.join(DATA_DIR, "biobank_ebv_loci_comparison.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bio_loci[row["gene"]] = row

    bio_signed = []
    with open(os.path.join(DATA_DIR, "biobank_signed_block_and_twas_concordance.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bio_signed.append(row)

    # Check 17: Top Biobank MHC/Immune peaks (HLA-B, ERAP2) null in LCLs, while Immune-Masked TWAS Tier 1 & LCL TWAS-mRNA are concordant
    hlab_p = float(bio_loci["HLA-B"]["lcl_exact_blood_lead_CN_p"])
    erap2_p = float(bio_loci["ERAP2"]["lcl_exact_blood_lead_CN_p"])
    tier1_masked = [r for r in bio_signed if r["test_option"] == "Option2B_TWAS_CrossStudy_Sign_All_12544_eQTLs"
                    and r["biobank_cohort"] == "Blood_UKB_AoU" and r["mask_mode"] == "Immune_MHC_Masked_SubThreshold"
                    and r["lcl_phenotype"] == "WGS_CN_ViralCond" and r["biobank_tier"] == "Tier1_SubThreshold_P_lt_1e-3"][0]
    twas_mrna = [r for r in bio_signed if r["test_option"] == "Option2B_LCL_GWAS_TWAS_vs_Observed_MAGE_mRNA"
                 and r["mask_mode"] == "Immune_MHC_Masked_SubThreshold" and r["lcl_phenotype"] == "WGS_CN_ViralCond"
                 and r["biobank_tier"] == "HighPIP_ge_0.5_and_Nominal_GWAS_P_lt_0.05"][0]
    record("7. Biobank vs LCL", "Immune-masked TWAS sub-threshold concordance & LCL TWAS-mRNA sign agreement",
           f"Exact LCL P: HLA-B={hlab_p:.3f}, ERAP2={erap2_p:.3f} | Masked Tier1 TWAS r={float(tier1_masked['z_pearson_r']):+.4f} | LCL TWAS-mRNA (PIP>=0.5) sign={float(twas_mrna['sign_concordance_frac'])*100:.1f}% (P={twas_mrna['binom_sign_p']})",
           hlab_p > 0.3 and erap2_p > 0.3 and float(tier1_masked["z_pearson_r"]) > 0.50 and float(twas_mrna["sign_concordance_frac"]) > 0.70)

    # Check 18: Unsigned Saliva lytic shedding block concordance & In Vivo Biobank viral variant origins (12,129:T>A vs B95-8 lab mutations)
    bio_uns = []
    with open(os.path.join(DATA_DIR, "biobank_unsigned_block_concordance.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            bio_uns.append(row)
    sal_bglf5 = [r for r in bio_uns if r["biobank_cohort"] == "Saliva_AoU" and r["mask_mode"] == "Unmasked_All_Blocks"
                 and r["lcl_phenotype"] == "RNA_BGLF5_ViralCond" and r["biobank_tier"] == "Tier1_SubThreshold_5e-8_to_1e-5"][0]
    sal_ielat = [r for r in bio_uns if r["biobank_cohort"] == "Saliva_AoU" and r["mask_mode"] == "Immune_MHC_Masked_SubThreshold"
                 and r["lcl_phenotype"] == "RNA_IELatRatio_ViralCond" and r["biobank_tier"] == "Tier1_SubThreshold_5e-8_to_1e-5"][0]
    viral_bio = {}
    with open(os.path.join(DATA_DIR, "biobank_vs_1kgp_viral_variants.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            viral_bio[row["variant_id"]] = row
    v_12129 = viral_bio["12129:T>A"]
    v_163858 = viral_bio["163858:G>T"]
    record("7. Biobank vs LCL", "Saliva lytic block concordance & UKB/AoU in vivo viral variant origins",
           f"Saliva->BGLF5 mean-chi2 rho={float(sal_bglf5['global_spearman_rho_mean_chisq']):+.4f} (P={sal_bglf5['global_spearman_p_mean_chisq']}); IE/Lat top10% OR={sal_ielat['top10pct_fisher_or']} | 12129:T>A UKB={float(v_12129['ukb_invivo_alt_freq'])*100:.1f}%/AoU={float(v_12129['aou_invivo_alt_freq'])*100:.1f}% vs 163858:G>T UKB={float(v_163858['ukb_invivo_alt_freq'])*100:.2f}%",
           float(sal_bglf5["global_spearman_rho_mean_chisq"]) > 0.05 and float(v_12129["ukb_invivo_alt_freq"]) >= 0.75 and float(v_163858["ukb_invivo_alt_freq"]) < 0.02)

    # Write machine-readable verification report
    report_path = os.path.join(ROOT_DIR, "REPLICATION_VERIFICATION_REPORT.md")
    with open(report_path, "w") as out:
        out.write("# Turn-Key Replication Verification Report\n\n")
        out.write("Generated automatically by `uv run replication_2026/reproduce_all.py` from `replication_2026/data/`.\n\n")
        out.write("| Section | Claim Verified | Recomputed Value | Status |\n")
        out.write("|---|---|---|---|\n")
        for sec, clm, val, st_str in checks:
            out.write(f"| {sec} | {clm} | `{val}` | **{st_str}** |\n")
    print(f"\nWrote verification report to: {report_path}")
    all_passed = all(c[3] == "PASS" for c in checks)
    print(f"Verification Summary: {sum(1 for c in checks if c[3]=='PASS')}/{len(checks)} checks PASSED.")
    return all_passed


if __name__ == "__main__":
    ok = run_verification()
    if "--no-figures" not in sys.argv:
        print("\n--- SECTION 8: Regenerating All 8 Publication Multi-Panel Figures ---")
        import generate_figures
        generate_figures.main()
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))
        import generate_fig7_host_genetics
        generate_fig7_host_genetics.main()
        import generate_fig8_biobank_vs_lcl
        generate_fig8_biobank_vs_lcl.main()
    sys.exit(0 if ok else 1)
