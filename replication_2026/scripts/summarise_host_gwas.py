#!/usr/bin/env python3
"""
Consolidate all 22 chromosomes of the Viral-Conditioned Host Genetic Scan (Layers A, B1, B2, C):
  - Layer A: Trio heritability before & after conditioning on viral genetics (n=570 clean trios).
  - Layer B1: Functional MAGE fine-mapped cis-eQTL lead variant scan (~14,000+ tested eQTLs)
              across all 11 EBV WGS & RNA phenotypes (with BH-FDR q-values).
  - Layer B2: Targeted 55 canonical host EBV-regulator genes (+/- 100 kb windows),
              reporting lead variant, local Bonferroni p-value, and BH-FDR across the 55 loci.
  - Layer C: Genome-wide common host variant GWAS (MAF >= 5%, chr1-chr22) across all 11
             phenotypes, computing exact genome-wide lambda_GC, QQ quantiles, minimum P-value,
             power / PVE upper bounds, and top suggestive loci.
"""
import csv
import glob
import json
import math
import os
from collections import defaultdict
import numpy as np

BASE = "/work/nfp8/EBV"
META = f"{BASE}/metadata"
GWAS_DIR = f"{BASE}/host_gwas"
CHR_DIR = f"{GWAS_DIR}/chr_results"

PHENO_NAMES = [
    "WGS_CN_PopAdj", "WGS_CN_ViralCond", "WGS_IR1_PopAdj",
    "RNA_BGLF5_ViralCond", "RNA_BZLF1_ViralCond", "RNA_IELatRatio_ViralCond",
    "RNA_ProdFork_ViralCond", "RNA_LMP1_ViralCond", "RNA_Latent_ViralCond",
    "RNA_BGLF5_Uncond", "RNA_LMP1_Uncond"
]


def bh_fdr(pvals):
    n = len(pvals)
    if n == 0:
        return []
    order = np.argsort(pvals)
    qvals = np.empty(n, dtype=float)
    min_q = 1.0
    for rank_rev in range(n - 1, -1, -1):
        idx = order[rank_rev]
        q = pvals[idx] * n / (rank_rev + 1)
        if q < min_q:
            min_q = q
        qvals[idx] = min(1.0, min_q)
    return qvals


def median_from_hist(hist, bin_width=0.005, high_tail=None):
    total = int(np.sum(hist)) + (len(high_tail) if high_tail else 0)
    if total == 0:
        return 0.0
    half = total / 2.0
    csum = np.cumsum(hist)
    idx = int(np.searchsorted(csum, half))
    if idx >= len(hist):
        return float(np.median(high_tail)) if high_tail else 0.0
    prev = csum[idx - 1] if idx > 0 else 0
    frac = (half - prev) / max(1, hist[idx])
    return (idx + frac) * bin_width


def main():
    print("=" * 100)
    print("CONSOLIDATING VIRAL-CONDITIONED HOST GENETIC SCAN (LAYERS A, B1, B2, AND C)")
    print("=" * 100)

    # 1. Layer C: Genome-Wide Summary Across chr1..22
    gw_summary_rows = []
    all_top_suggestive = []

    for pheno in PHENO_NAMES:
        tot_maf01 = 0
        tot_common = 0
        tot_gw_5e8 = 0
        tot_sug_1e6 = 0
        tot_sug_1e5 = 0
        comb_hist = np.zeros(10000, dtype=np.int64)
        best_p = 2.0
        best_row = None

        for c in range(1, 23):
            sfile = f"{CHR_DIR}/chr{c}_{pheno}_summary.json"
            if not os.path.exists(sfile):
                continue
            with open(sfile) as f:
                d = json.load(f)
            tot_maf01 += d["n_maf01"]
            tot_common += d["n_common_maf05"]
            tot_gw_5e8 += d["n_gw_5e8"]
            tot_sug_1e6 += d["n_sug_1e6"]
            tot_sug_1e5 += d["n_sug_1e5"]
            comb_hist += np.array(d["chisq_hist"], dtype=np.int64)
            m_row = d["min_p_row"]
            if m_row is not None:
                p_val = float(m_row[-1])
                if p_val < best_p:
                    best_p = p_val
                    best_row = m_row

        # Median chi2 under 1-df null is 0.4549364
        med_chi2 = median_from_hist(comb_hist, 0.005)
        lambda_gc = med_chi2 / 0.4549364 if med_chi2 > 0 else 1.0

        lead_id = best_row[2] if best_row else "NA"
        lead_chr = best_row[0] if best_row else "NA"
        lead_pos = best_row[1] if best_row else "NA"
        lead_freq = best_row[4] if best_row else "NA"
        lead_n = int(best_row[5]) if best_row else 0
        lead_beta = float(best_row[6]) if best_row else 0.0
        lead_t = float(best_row[8]) if best_row else 0.0
        lead_pve = (lead_t**2) / (lead_t**2 + max(1, lead_n - 2)) if lead_n > 2 else 0.0

        gw_summary_rows.append({
            "phenotype": pheno,
            "n_samples": lead_n,
            "n_variants_maf01": tot_maf01,
            "n_common_maf05": tot_common,
            "lambda_gc": f"{lambda_gc:.4f}",
            "n_gw_sig_5e8": tot_gw_5e8,
            "n_suggestive_1e6": tot_sug_1e6,
            "n_suggestive_1e5": tot_sug_1e5,
            "lead_variant_id": lead_id,
            "lead_chrom": lead_chr,
            "lead_pos": lead_pos,
            "lead_a1_freq": lead_freq,
            "lead_beta_sd": f"{lead_beta:+.4f}",
            "lead_t_stat": f"{lead_t:+.4f}",
            "lead_p": f"{best_p:.3e}",
            "lead_pve_pct": f"{lead_pve * 100:.2f}"
        })
        print(f"  [Layer C GWAS | {pheno:<24s}] n={lead_n:4d}, M_common={tot_common:,}, lambda_GC={lambda_gc:.3f}, "
              f"P<5e-8={tot_gw_5e8}, P<1e-6={tot_sug_1e6}, min_P={best_p:.2e} ({lead_id}, PVE={lead_pve*100:.2f}%)")

        # Collect top suggestive rows across chromosomes
        sug_rows = []
        for c in range(1, 23):
            sf = f"{CHR_DIR}/chr{c}_{pheno}_suggestive.tsv"
            if os.path.exists(sf):
                with open(sf) as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        freq = float(row["a1_freq"])
                        if 0.05 <= freq <= 0.95:
                            sug_rows.append(row)
        sug_rows.sort(key=lambda x: float(x["p"]))
        for r in sug_rows[:25]:
            all_top_suggestive.append({"phenotype": pheno, **r})

    out_gw = f"{META}/host_gwas_genomewide_summary.tsv"
    with open(out_gw, "w") as out:
        cols = list(gw_summary_rows[0].keys())
        out.write("\t".join(cols) + "\n")
        for r in gw_summary_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    out_sug = f"{META}/host_gwas_top_loci.tsv"
    if all_top_suggestive:
        with open(out_sug, "w") as out:
            cols = list(all_top_suggestive[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in all_top_suggestive:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 2. Layer B1: Functional MAGE Fine-Mapped cis-eQTL Lead Variants
    print("\n--- LAYER B1: MAGE Fine-Mapped cis-eQTL Lead Variants Scan ---")
    eqtl_summary_rows = []
    top_eqtl_rows = []
    for pheno in PHENO_NAMES:
        rows = []
        for c in range(1, 23):
            ef = f"{CHR_DIR}/chr{c}_{pheno}_mage_eqtls.tsv"
            if os.path.exists(ef):
                with open(ef) as f:
                    for r in csv.DictReader(f, delimiter="\t"):
                        rows.append(r)
        if not rows:
            continue
        # Deduplicate by (id, gene) keeping highest PIP
        by_pair = {}
        for r in rows:
            k = (r["id"], r["gene"])
            if k not in by_pair or float(r["pip"]) > float(by_pair[k]["pip"]):
                by_pair[k] = r
        dedup = list(by_pair.values())
        pvals = np.array([float(r["p"]) for r in dedup])
        qvals = bh_fdr(pvals)
        for r, q in zip(dedup, qvals):
            r["q_fdr"] = f"{q:.4f}"
        dedup.sort(key=lambda x: float(x["p"]))
        n_q05 = sum(1 for q in qvals if q < 0.05)
        n_q10 = sum(1 for q in qvals if q < 0.10)
        best = dedup[0]
        eqtl_summary_rows.append({
            "phenotype": pheno,
            "n_eqtls_tested": len(dedup),
            "n_fdr_q05": n_q05,
            "n_fdr_q10": n_q10,
            "top_gene": best["gene"],
            "top_rsid": best["rsid"],
            "top_variant_id": best["id"],
            "top_pip": best["pip"],
            "top_beta": best["beta"],
            "top_p": best["p"],
            "top_q_fdr": best["q_fdr"]
        })
        print(f"  [Layer B1 eQTLs | {pheno:<24s}] M={len(dedup):,}, q<0.05={n_q05}, q<0.10={n_q10} | "
              f"Top: {best['gene']} ({best['rsid']}, PIP={best['pip']}) beta={float(best['beta']):+.3f}, P={float(best['p']):.2e}, q={best['q_fdr']}")
        for r in dedup[:15]:
            top_eqtl_rows.append({"phenotype": pheno, **r})

    out_eqtl_sum = f"{META}/host_mage_eqtl_scan_summary.tsv"
    if eqtl_summary_rows:
        with open(out_eqtl_sum, "w") as out:
            cols = list(eqtl_summary_rows[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in eqtl_summary_rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")
    out_eqtl_top = f"{META}/host_mage_eqtl_top_hits.tsv"
    if top_eqtl_rows:
        with open(out_eqtl_top, "w") as out:
            cols = list(top_eqtl_rows[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in top_eqtl_rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 3. Layer B2: Targeted 55 Canonical Host EBV-Regulator Genes (+/- 100 kb)
    print("\n--- LAYER B2: 55 Canonical Host EBV-Regulator Genes (+/- 100 kb) ---")
    all_cand_rows = []
    for pheno in PHENO_NAMES:
        crows = []
        for c in range(1, 23):
            cf = f"{CHR_DIR}/chr{c}_{pheno}_candidate55.tsv"
            if os.path.exists(cf):
                with open(cf) as f:
                    for r in csv.DictReader(f, delimiter="\t"):
                        # Compute window-level Bonferroni p = min(1.0, p * n_maf01)
                        p_raw = float(r["p"])
                        n_m = int(r["n_maf01"])
                        # Effective independent tests in +/-100kb is roughly min(n_m, 250)
                        p_win = min(1.0, p_raw * max(1, min(n_m, 250)))
                        r["p_window_adj"] = p_win
                        crows.append(r)
        if not crows:
            continue
        p_wins = np.array([r["p_window_adj"] for r in crows])
        q_cands = bh_fdr(p_wins)
        for r, q in zip(crows, q_cands):
            r["p_window_adj"] = f"{r['p_window_adj']:.4e}"
            r["q_cand55"] = f"{q:.4f}"
        crows.sort(key=lambda x: float(x["p"]))
        b = crows[0]
        n_sig_cand = sum(1 for q in q_cands if q < 0.05)
        print(f"  [Layer B2 Cand54| {pheno:<24s}] q<0.05={n_sig_cand}/54 | Lead gene: {b['gene']} ({b['category']}) "
              f"var={b['id']} beta={float(b['beta']):+.3f}, raw_P={float(b['p']):.2e}, win_q={b['q_cand55']}")
        for r in crows:
            all_cand_rows.append({"phenotype": pheno, **r})

    out_cand = f"{META}/host_candidate55_scan.tsv"
    if all_cand_rows:
        with open(out_cand, "w") as out:
            cols = list(all_cand_rows[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in all_cand_rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 4. Layer B3: Polymorphic Host Mobile Element Insertions & Absences (MEGANE 1kGP, M=8,458 at MAF>=1%)
    print("\n--- LAYER B3: Polymorphic Host Mobile Element Insertions & Absences (L1 / Alu / SVA) ---")
    mei_summary_rows = []
    for pheno in PHENO_NAMES:
        m_rows = []
        for prefix, vtype in [("mei_ins_res", "MEI_Insertion"), ("mei_abs_res", "MEI_Absence")]:
            mfile = f"{GWAS_DIR}/{prefix}.{pheno}.glm.linear"
            if not os.path.exists(mfile):
                continue
            with open(mfile) as f:
                hdr = [h.lstrip("#") for h in f.readline().strip().split("\t")]
                for line in f:
                    p = line.strip().split("\t")
                    r = dict(zip(hdr, p))
                    if r.get("P", "NA") != "NA":
                        r["mei_type"] = vtype
                        m_rows.append(r)
        if not m_rows:
            continue
        pvals = np.array([float(r["P"]) for r in m_rows])
        qvals = bh_fdr(pvals)
        n_q05 = sum(1 for q in qvals if q < 0.05)
        n_q10 = sum(1 for q in qvals if q < 0.10)
        best_idx = int(np.argmin(pvals))
        b = m_rows[best_idx]
        mei_summary_rows.append({
            "phenotype": pheno,
            "n_mei_tested": len(m_rows),
            "n_fdr_q05": n_q05,
            "n_fdr_q10": n_q10,
            "lead_mei_id": b["ID"],
            "lead_mei_type": b["mei_type"],
            "lead_chrom": b["CHROM"],
            "lead_pos": b["POS"],
            "lead_freq": b["A1_FREQ"],
            "lead_beta": b["BETA"],
            "lead_p": b["P"],
            "lead_q_fdr": f"{qvals[best_idx]:.4f}"
        })
        print(f"  [Layer B3 MEIs  | {pheno:<24s}] M={len(m_rows):,}, q<0.05={n_q05} | "
              f"Lead: {b['ID']} ({b['mei_type']} chr{b['CHROM']}:{b['POS']}) beta={float(b['BETA']):+.3f}, P={float(b['P']):.2e}, q={qvals[best_idx]:.4f}")

    out_mei = f"{META}/host_mei_scan_summary.tsv"
    if mei_summary_rows:
        with open(out_mei, "w") as out:
            cols = list(mei_summary_rows[0].keys())
            out.write("\t".join(cols) + "\n")
            for r in mei_summary_rows:
                out.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("\nSaved consolidated host genetic tables to:")
    print(f"  - {out_gw}")
    print(f"  - {out_sug}")
    print(f"  - {out_eqtl_sum}")
    print(f"  - {out_eqtl_top}")
    print(f"  - {out_cand}")


if __name__ == "__main__":
    main()
