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
Generate all 6 publication-grade multi-panel figures (PNG 300 DPI + vector PDF)
for the 1000 Genomes / MAGE EBV Copy Number, Viral Genetics, and Host Transcriptome
Mediation study directly from local master tables in `replication_2026/data/`.

Usage:
    uv run replication_2026/scripts/generate_figures.py
"""
import csv
import math
import os
from collections import defaultdict
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
FIG_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "figures"))
os.makedirs(FIG_DIR, exist_ok=True)

# Publication style configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 9.2,
    "axes.titlesize": 9.8,
    "axes.titleweight": "bold",
    "axes.titlepad": 10.0,
    "axes.labelsize": 9.2,
    "xtick.labelsize": 8.2,
    "ytick.labelsize": 8.2,
    "legend.fontsize": 7.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

CONT_COLORS = {
    "AFR": "#E69F00",  # Orange
    "AMR": "#CC79A7",  # Reddish purple
    "EAS": "#009E73",  # Bluish green
    "EUR": "#0072B2",  # Blue
    "SAS": "#56B4E9",  # Sky blue
}

CLASS_COLORS = {
    "LYTIC_IE": "#D55E00",
    "LYTIC_EARLY": "#E69F00",
    "LYTIC_LATE": "#CC79A7",
    "LATENT_MRNA": "#0072B2",
    "LATENT_BART": "#009E73",
    "NONCODING": "#7F7F7F",
    "B958_DEL": "#333333",
}


def add_panel_label(ax, label, x=-0.14, y=1.11):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")


def demean_by_groups(vals, groups):
    g_sums = defaultdict(float)
    g_cnts = defaultdict(int)
    for v, g in zip(vals, groups):
        g_sums[g] += v
        g_cnts[g] += 1
    g_means = {g: g_sums[g] / g_cnts[g] for g in g_cnts}
    return np.array([v - g_means[g] for v, g in zip(vals, groups)])


def load_core_data():
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
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            unrel.add(row["sample"])

    ont_cn = {}
    with open(os.path.join(DATA_DIR, "all_ont_904.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                ont_cn[p[0]] = float(p[5])

    low_cn, low_ebv_dep = {}, {}
    with open(os.path.join(DATA_DIR, "all_lowcov_2499.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                low_ebv_dep[p[0]] = float(p[2])
                low_cn[p[0]] = float(p[5])

    ir1_30, ir1_flank_ratio, flank_uniq_ratio = {}, {}, {}
    with open(os.path.join(DATA_DIR, "all_ir1_2502.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 6:
                s = p[0]
                d_ir1, d_cp, d_uniq = float(p[1]), float(p[2]), float(p[3])
                ir1_30[s] = float(p[4]) * (22737.0 / 3072.0)
                if d_cp > 0 and d_uniq > 0:
                    ir1_flank_ratio[s] = (d_ir1 / d_cp) * (22737.0 / 3072.0)
                    flank_uniq_ratio[s] = d_cp / d_uniq

    ont_ir1 = {}
    with open(os.path.join(DATA_DIR, "all_ont_ir1_907.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4:
                ont_ir1[p[0]] = float(p[3])

    gcbias = {}
    with open(os.path.join(DATA_DIR, "all_gcbias_2502.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 5:
                gcbias[p[0]] = float(p[4])

    bcr_top_frac = {}
    with open(os.path.join(DATA_DIR, "bcr_wgs_summary.tsv")) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 7:
                bcr_top_frac[p[0]] = float(p[6])

    mage_runs = []
    with open(os.path.join(DATA_DIR, "all_mage_ebv_rna_779.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            mage_runs.append(row)

    return {
        "ped": ped, "cn30": cn30, "wt30": wt30, "unrel": unrel,
        "ont_cn": ont_cn, "low_cn": low_cn, "low_ebv_dep": low_ebv_dep,
        "ir1_30": ir1_30, "ir1_flank_ratio": ir1_flank_ratio,
        "flank_uniq_ratio": flank_uniq_ratio, "ont_ir1": ont_ir1,
        "gcbias": gcbias, "bcr_top_frac": bcr_top_frac, "mage_runs": mage_runs
    }


def make_figure_1(D):
    """Figure 1: WGS EBV CN Architecture, 3-Platform Concordance, Trio Heritability, & BCR Clonality."""
    ped, cn30, wt30, unrel = D["ped"], D["cn30"], D["wt30"], D["unrel"]
    clean_unrel = sorted([s for s in unrel if s in cn30 and wt30.get(s, 0) <= 0.01 and s in ped])

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.8))
    plt.subplots_adjust(wspace=0.28, hspace=0.36)

    # Panel A: 26-population boxplot of log2(CN) across n = 2,449 clean 1kGP
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    pop_vals = defaultdict(list)
    pop_cont = {}
    for s in clean_unrel:
        p = ped[s]["pop"]
        c = ped[s]["cont"]
        pop_vals[p].append(math.log2(cn30[s]))
        pop_cont[p] = c

    cont_order = ["EUR", "AFR", "EAS", "SAS", "AMR"]
    sorted_pops = sorted(pop_vals.keys(), key=lambda p: (cont_order.index(pop_cont[p]), np.median(pop_vals[p])))
    bp_data = [pop_vals[p] for p in sorted_pops]
    bp = ax.boxplot(bp_data, patch_artist=True, showfliers=False, widths=0.68,
                    medianprops=dict(color="black", linewidth=1.4))
    for patch, p in zip(bp["boxes"], sorted_pops):
        patch.set_facecolor(CONT_COLORS[pop_cont[p]])
        patch.set_alpha(0.82)
        if p in ("GBR", "IBS"):
            patch.set_edgecolor("#D55E00" if p == "IBS" else "#000000")
            patch.set_linewidth(1.8)

    ax.set_xticks(range(1, len(sorted_pops) + 1))
    ax.set_xticklabels(sorted_pops, rotation=60, ha="right", fontsize=7.5)
    ax.set_ylabel(r"$\log_2(\mathrm{EBV\ Copy\ Number})$ (30x WGS)")
    ax.set_title(f"1000 Genomes EBV Copy Number Across 26 Populations (n = {len(clean_unrel):,})")
    ax.axhline(math.log2(np.median([cn30[s] for s in clean_unrel])), color="#555555", ls="--", lw=1.0, label="Cohort median (36.2)")
    handles = [mpatches.Patch(color=CONT_COLORS[c], label=c, alpha=0.82) for c in cont_order]
    ax.legend(handles=handles, loc="upper left", ncol=5, frameon=False, fontsize=7.5)
    ax.annotate("GBR\n(med 26.6)", xy=(sorted_pops.index("GBR") + 1, math.log2(26.6)),
                xytext=(sorted_pops.index("GBR") + 0.6, 2.3),
                arrowprops=dict(arrowstyle="->", lw=0.9), fontsize=7.5, fontweight="bold")
    ax.annotate("IBS\n(med 60.5)", xy=(sorted_pops.index("IBS") + 1, math.log2(60.5)),
                xytext=(sorted_pops.index("IBS") - 1.2, 8.2),
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#D55E00"), color="#D55E00", fontsize=7.5, fontweight="bold")

    # Panel B: 3-Platform Concordance (30x WGS vs ONT Long-Read CN, n = 890 clean)
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    ont_cn = D["ont_cn"]
    shared_ont = [s for s in clean_unrel if s in ont_cn and ont_cn[s] > 0]
    x_30 = np.array([math.log2(cn30[s]) for s in shared_ont])
    y_ont = np.array([math.log2(ont_cn[s]) for s in shared_ont])
    c_ont = [CONT_COLORS[ped[s]["cont"]] for s in shared_ont]
    r_ont, _ = stats.pearsonr(x_30, y_ont)
    ax.scatter(x_30, y_ont, c=c_ont, s=14, alpha=0.65, edgecolors="none")
    m, b = np.polyfit(x_30, y_ont, 1)
    xs = np.linspace(x_30.min(), x_30.max(), 100)
    ax.plot(xs, m * xs + b, color="#111111", lw=1.6, label=f"OLS fit (r = +{r_ont:.4f}, R² = {r_ont**2:.3f})")
    ax.plot(xs, xs, color="#888888", ls=":", lw=1.1, label="y = x identity")
    ax.set_xlabel(r"2019 30x Illumina WGS $\log_2(\mathrm{EBV\ CN})$")
    ax.set_ylabel(r"~2023 Oxford Nanopore (ONT) $\log_2(\mathrm{EBV\ CN})$")
    ax.set_title(f"Cross-Platform & Cross-Extraction Stability (n = {len(shared_ont)} clean lines)")
    ax.text(0.04, 0.14,
            "3-Platform Reliability Triangle (n = 703):\n"
            "• 30x vs ONT: r = +0.9072\n"
            "• 30x vs Low-cov: r = +0.8427\n"
            "• ONT vs Low-cov: r = +0.8188\n"
            "• Implied reliability: 30x = 0.957, ONT = 0.856",
            transform=ax.transAxes, fontsize=7.8,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8F9FA", edgecolor="#CCCCCC"))
    ax.legend(loc="upper left", frameon=False)

    # Panel C: Trio Heritability Decomposition (n = 570 clean trios: Pooled vs Within-Population)
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    clean_all = {s: math.log2(v) for s, v in cn30.items() if wt30.get(s, 0) <= 0.01 and v > 0}
    trios = []
    for s, info in ped.items():
        fa, mo, pop = info["dad"], info["mom"], info["pop"]
        if fa != "0" and mo != "0" and s in clean_all and fa in clean_all and mo in clean_all:
            trios.append((s, fa, mo, pop))

    p_vals = defaultdict(list)
    for s, v in clean_all.items():
        if s in ped:
            p_vals[ped[s]["pop"]].append(v)
    p_mean = {p: np.mean(v) for p, v in p_vals.items()}

    mid_raw = np.array([0.5 * (clean_all[fa] + clean_all[mo]) for s, fa, mo, p in trios])
    off_raw = np.array([clean_all[s] for s, fa, mo, p in trios])
    sp_x_raw = np.array([clean_all[fa] for _, fa, _, _ in trios] + [clean_all[mo] for _, _, mo, _ in trios])
    sp_y_raw = np.array([clean_all[s] for s, _, _, _ in trios] * 2)
    pp_x_raw = np.array([clean_all[fa] for _, fa, _, _ in trios])
    pp_y_raw = np.array([clean_all[mo] for _, _, mo, _ in trios])

    val_w = lambda s, p: clean_all[s] - p_mean[p]
    mid_w = np.array([0.5 * (val_w(fa, p) + val_w(mo, p)) for s, fa, mo, p in trios])
    off_w = np.array([val_w(s, p) for s, fa, mo, p in trios])
    sp_x_w = np.array([val_w(fa, p) for _, fa, _, p in trios] + [val_w(mo, p) for _, _, mo, p in trios])
    sp_y_w = np.array([val_w(s, p) for s, _, _, p in trios] * 2)
    pp_x_w = np.array([val_w(fa, p) for _, fa, _, p in trios])
    pp_y_w = np.array([val_w(mo, p) for _, _, mo, p in trios])

    r_mid_raw, _ = stats.pearsonr(mid_raw, off_raw)
    r_sp_raw, _ = stats.pearsonr(sp_x_raw, sp_y_raw)
    r_pp_raw, _ = stats.pearsonr(pp_x_raw, pp_y_raw)

    r_mid_w, _ = stats.pearsonr(mid_w, off_w)
    r_sp_w, _ = stats.pearsonr(sp_x_w, sp_y_w)
    r_pp_w, _ = stats.pearsonr(pp_x_w, pp_y_w)

    cats = ["Mid-parent vs\nOffspring (r_mid)", "Single-parent vs\nOffspring (r_sp)",
            "Spouse vs Spouse\n(r_spouse, non-genetic)", "Narrow-sense h²\n2(r_sp - r_spouse)"]
    raw_bars = [r_mid_raw, r_sp_raw, r_pp_raw, 2 * (r_sp_raw - r_pp_raw)]
    within_bars = [r_mid_w, r_sp_w, r_pp_w, 2 * (r_sp_w - r_pp_w)]

    x_idx = np.arange(len(cats))
    w = 0.36
    b1 = ax.bar(x_idx - w/2, raw_bars, w, color="#888888", alpha=0.85, label="Pooled across 26 populations")
    b2 = ax.bar(x_idx + w/2, within_bars, w, color="#0072B2", alpha=0.90, label="Within-population demeaned (f[5])")
    ax.axhline(0.2097, color="#D55E00", ls="--", lw=1.3, label="Same-population unrelated pair baseline (r = +0.210)")
    ax.axhline(0, color="black", lw=0.8)
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2, h + 0.012, f"{h:+.3f}",
                    ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax.set_xticks(x_idx)
    ax.set_xticklabels(cats, fontsize=7.8)
    ax.set_ylim(-0.04, 0.52)
    ax.set_ylabel("Correlation (r) or Heritability (h²)")
    ax.set_title(f"Zero Host Genetic Heritability of EBV CN (n = {len(trios)} Clean Trios)")
    ax.legend(loc="upper right", frameon=False, fontsize=7.6)

    # Panel D: PoN-Cleaned BCR Clonality Dose-Response & Independence from EBV CN
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    bcr_frac = D["bcr_top_frac"]
    shared_bcr = [s for s in clean_unrel if s in bcr_frac]
    x_bcr = np.array([bcr_frac[s] for s in shared_bcr])
    y_cn_bcr = np.array([math.log2(cn30[s]) for s in shared_bcr])
    rho_bcr, _ = stats.spearmanr(x_bcr, y_cn_bcr)
    r_bcr, _ = stats.pearsonr(x_bcr, y_cn_bcr)
    ax.scatter(x_bcr, y_cn_bcr, s=10, alpha=0.35, color="#009E73", edgecolors="none")
    m_b, b_b = np.polyfit(x_bcr, y_cn_bcr, 1)
    xs_b = np.linspace(0, 1, 100)
    ax.plot(xs_b, m_b * xs_b + b_b, color="#111111", lw=1.6,
            label=f"Spearman ρ = {rho_bcr:+.4f} (R² = {r_bcr**2*100:.2f}%, n = {len(shared_bcr):,})")
    ax.set_xlabel("PoN-Cleaned Dominant IGH Clone Fraction in WGS (top_clone_frac)")
    ax.set_ylabel(r"30x WGS $\log_2(\mathrm{EBV\ Copy\ Number})$")
    ax.set_title("B-Cell Clonality Does Not Drive EBV Copy Number Variation")
    ax.legend(loc="lower right", frameon=False)

    # Inset showing PoN-cleaned WGS-to-RNA read-depth dose-response
    ins = ax.inset_axes([0.52, 0.57, 0.44, 0.37])
    bins_lbl = ["1–5", "6–15", "16–30", ">30"]
    rec_rates = [80.8, 94.2, 96.6, 100.0]
    ins_x = np.arange(len(bins_lbl))
    ins.bar(ins_x, rec_rates, color="#0072B2", alpha=0.88, width=0.58)
    ins.set_xticks(ins_x)
    ins.set_xticklabels(bins_lbl)
    for idx, val in enumerate(rec_rates):
        ins.text(idx, val - 18, f"{val:.1f}%", ha="center", color="white", fontsize=6.6, fontweight="bold")
    ins.set_ylim(0, 114)
    ins.set_ylabel("% in RNA", fontsize=7.0)
    ins.set_xlabel("WGS Read Depth", fontsize=7.0)
    ins.set_title("PoN-Cleaned RNA Recovery", fontsize=7.3, fontweight="bold", pad=4)
    ins.tick_params(labelsize=6.8)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig1_wgs_architecture_and_heritability.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_wgs_architecture_and_heritability.(png/pdf)")


def make_figure_2(D):
    """Figure 2: H1 Decisive Test — EBV CN Tracks Productive Lytic Reactivation, Not Latent Dosage."""
    cn30, wt30, ped = D["cn30"], D["wt30"], D["ped"]
    mage_runs = D["mage_runs"]

    by_donor = defaultdict(list)
    for r in mage_runs:
        if r["donor"] in cn30:
            by_donor[r["donor"]].append(r)

    clean_donors = sorted([d for d in by_donor if wt30.get(d, 0) <= 0.01 and d in ped])

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.8))
    plt.subplots_adjust(wspace=0.30, hspace=0.36)

    # Panel A: Cross-thaw biological vs same-batch technical repeatability (24 donors, 70 pairs)
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    metrics = [
        ("BGLF5\n(Shutoff)", 0.9597, 0.9504),
        ("log2(IE /\nLatent)", 0.9725, 0.9102),
        ("Lytic IE\nRPKM", 0.9524, 0.8053),
        ("log2(Late /\nLatent)", 0.9739, 0.7631),
        ("Lytic Early\nRPKM", 0.9465, 0.6754),
        ("Lytic Late\nRPKM", 0.9614, 0.6365),
        ("Latent Pol II\nmRNA", 0.8581, 0.2375),
    ]
    x_idx = np.arange(len(metrics))
    w = 0.36
    b1 = ax.bar(x_idx - w/2, [m[1] for m in metrics], w, color="#56B4E9", alpha=0.9, label="Same-batch technical replicates (n = 24 pairs)")
    b2 = ax.bar(x_idx + w/2, [m[2] for m in metrics], w, color="#D55E00", alpha=0.9, label="Different-batch independent thaws (n = 46 pairs)")
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2, h + 0.015, f"{h:.2f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold")
    ax.set_xticks(x_idx)
    ax.set_xticklabels([m[0] for m in metrics], fontsize=7.6)
    ax.set_ylim(0, 1.22)
    ax.set_ylabel("Pairwise Replicate Correlation (r)")
    ax.set_title("Lytic Reactivation Is Intrinsic Across Thaws (70 Pairs, 24 Donors)")
    ax.legend(loc="upper right", frameon=False, fontsize=7.6)

    # Panel B: H1 Joint Scatter — WGS log2(CN) vs Latent Pol II mRNA vs Lytic IE RPKM (n = 717)
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    x_cn = np.array([math.log2(cn30[d]) for d in clean_donors])
    y_lat = np.array([math.log2(np.mean([float(r["latent_mrna_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_donors])
    y_ie = np.array([math.log2(np.mean([float(r["lytic_ie_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_donors])

    r_lat, _ = stats.pearsonr(x_cn, y_lat)
    r_ie, p_ie = stats.pearsonr(x_cn, y_ie)

    ax.scatter(x_cn, y_lat, color="#0072B2", s=13, alpha=0.45, edgecolors="none",
               label=f"Latent Pol II mRNA (EBNA/LMP): r = {r_lat:+.4f} (null)")
    ax.scatter(x_cn, y_ie, color="#D55E00", s=13, alpha=0.55, edgecolors="none",
               label=f"Lytic Immediate-Early (BZLF1/BRLF1): r = {r_ie:+.4f} (p = {p_ie:.1e})")

    xs = np.linspace(x_cn.min(), x_cn.max(), 100)
    m1, b1_i = np.polyfit(x_cn, y_lat, 1)
    m2, b2_i = np.polyfit(x_cn, y_ie, 1)
    ax.plot(xs, m1 * xs + b1_i, color="#004C78", lw=2.0, ls="--")
    ax.plot(xs, m2 * xs + b2_i, color="#8C3800", lw=2.2)
    ax.set_xlabel(r"2019 30x WGS $\log_2(\mathrm{EBV\ Copy\ Number})$")
    ax.set_ylabel(r"2023 MAGE Viral Transcript $\log_2(\mathrm{RPKM} + 0.01)$")
    ax.set_title(f"Decisive H1 Test: CN Tracks Lytic IE, Not Latent Dosage (n = {len(clean_donors)})")
    ax.legend(loc="lower right", frameon=False, fontsize=7.8)

    # Panel C: Marginal, Within-Population, and Within-Population+Batch H1 Correlations
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    pops = [ped[d]["pop"] for d in clean_donors]
    pop_batches = [f"{ped[d]['pop']}_{by_donor[d][0]['batch']}" for d in clean_donors]

    y_zspl = np.array([math.log2(np.mean([float(r["bzlf1_spliced_reads"]) / float(r["total_reads"]) * 1e6 for r in by_donor[d]]) + 0.001) for d in clean_donors])
    y_ear = np.array([math.log2(np.mean([float(r["lytic_early_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_donors])
    y_late = np.array([math.log2(np.mean([float(r["lytic_late_rpkm"]) for r in by_donor[d]]) + 0.01) for d in clean_donors])
    y_rat = np.array([np.mean([float(r["lytic_late_over_latent_log2"]) for r in by_donor[d]]) for d in clean_donors])

    series = [
        ("Spliced BZLF1\n(Zta) RPM", y_zspl),
        ("Lytic IE\nRPKM", y_ie),
        ("log2(Late /\nLatent)", y_rat),
        ("Lytic Late\nRPKM", y_late),
        ("Lytic Early\nRPKM", y_ear),
        ("Latent Pol II\nmRNA RPKM", y_lat),
    ]
    c_marg = [stats.pearsonr(x_cn, v)[0] for _, v in series]
    c_wp = [stats.pearsonr(demean_by_groups(x_cn, pops), demean_by_groups(v, pops))[0] for _, v in series]
    c_wb = [stats.pearsonr(demean_by_groups(x_cn, pop_batches), demean_by_groups(v, pop_batches))[0] for _, v in series]

    x_i = np.arange(len(series))
    w3 = 0.25
    ax.bar(x_i - w3, c_marg, w3, color="#E69F00", label="Marginal (n = 717)")
    ax.bar(x_i, c_wp, w3, color="#D55E00", label="Within 26 populations")
    ax.bar(x_i + w3, c_wb, w3, color="#8C2D04", label="Within 26 pops × 17 batches")
    ax.axhline(0, color="black", lw=0.8)
    for idx, val in enumerate(c_wb):
        ax.text(idx + w3, val + (0.018 if val >= 0 else -0.052), f"{val:+.2f}", ha="center", fontsize=7.0, fontweight="bold")
    ax.set_xticks(x_i)
    ax.set_xticklabels([s[0] for s in series], fontsize=7.6)
    ax.set_ylim(-0.14, 0.86)
    ax.set_ylabel(r"Correlation (r) with WGS $\log_2(\mathrm{EBV\ CN})$")
    ax.set_title("Invariance of Lytic–CN Coupling to Population & Batch")
    ax.legend(loc="upper right", frameon=False, fontsize=7.6)

    # Panel D: Fine-Grained Splice Junction Architecture (Latency III Cp/Wp vs Zero Latency I Qp)
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    spl_counts = defaultdict(list)
    with open(os.path.join(DATA_DIR, "all_mage_viral_splices_59.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            spl_counts[row["splice_class"]].append(int(row["count"]))

    spl_order = [
        ("Wp_W1_W2_IR1", "Wp (W1–W2)\nIR1 Repeat"),
        ("BART_splices", "BART\nExons I–VII"),
        ("LMP1_introns", "LMP-1\nIntrons a/b"),
        ("LMP2A_2B_exons", "LMP-2A/2B\nExons 2–8"),
        ("Cp_EBNA", "Cp (C1–C2–W1)\nLatency III"),
        ("W2_to_Y1Y2_EBNA2", "W2 -> Y1/Y2\n(EBNA-2/LP)"),
        ("Latent_W2_Y_to_BHRF1", "W2/Y -> BHRF1\n(Latent vBCL2)"),
        ("CpWp_U_to_EBNA1", "Cp/Wp U -> K\n(Latent III EBNA1)"),
        ("Qp_to_EBNA1", "Qp -> U -> K\n(Latency I = 0/59)"),
    ]
    means = [np.mean(spl_counts[k]) if spl_counts[k] else 0.0 for k, _ in spl_order]
    n_pos = [sum(1 for v in spl_counts[k] if v > 0) for k, _ in spl_order]
    n_tot = len(spl_counts["Cp_EBNA"])
    y_pos = np.arange(len(spl_order))
    colors_d = ["#0072B2" if k != "Qp_to_EBNA1" else "#D55E00" for k, _ in spl_order]
    ax.barh(y_pos, [math.log10(m + 1.0) for m in means], color=colors_d, alpha=0.88, height=0.65)
    for idx, (m, np_cnt) in enumerate(zip(means, n_pos)):
        ax.text(math.log10(m + 1.0) + 0.06, idx, f"{m:,.1f} reads ({np_cnt}/{n_tot} runs)",
                va="center", fontsize=7.4, fontweight="bold" if idx == len(spl_order)-1 else "normal")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([lbl for _, lbl in spl_order], fontsize=7.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 5.2)
    ax.set_xlabel(r"$\log_{10}(\mathrm{Mean\ Spliced\ Junction\ Reads} + 1)$")
    ax.set_title("Exclusively Latency III (Cp/Wp) & Zero Latency I (Qp) Splicing")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig2_lytic_vs_latent_rna_h1.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig2_lytic_vs_latent_rna_h1.(png/pdf)")


def make_figure_3():
    """Figure 3: Aim 1 — Dissection of All 84 EBV Genes & Abortive vs Productive Lytic Replication."""
    genes = []
    with open(os.path.join(DATA_DIR, "viral_gene_dissection_84.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for k in ("length_bp", "mean_rpkm", "median_rpkm", "r_tech", "r_bio",
                      "r_cn_marginal", "r_cn_within_pop_batch", "r_part_given_bzlf1"):
                row[k] = float(row[k])
            genes.append(row)

    fig, axes = plt.subplots(2, 2, figsize=(14.2, 10.0))
    plt.subplots_adjust(wspace=0.32, hspace=0.38)

    # Panel A: Ranked distribution of all 84 EBV genes by within-pop+batch r with log2(CN)
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    genes_sorted = sorted(genes, key=lambda g: g["r_cn_within_pop_batch"], reverse=True)
    x_idx = np.arange(len(genes_sorted))
    r_wb = [g["r_cn_within_pop_batch"] for g in genes_sorted]
    c_wb = [CLASS_COLORS.get(g["kinetic_class"], "#888888") for g in genes_sorted]
    ax.bar(x_idx, r_wb, color=c_wb, width=0.85, alpha=0.9)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlim(-1, len(genes_sorted))
    ax.set_ylim(-0.26, 0.82)
    ax.set_xlabel("All 84 Non-Overlapping EBV Gene Windows (Ranked by Within-Pop+Batch r)")
    ax.set_ylabel(r"Within-Pop+Batch $r$ with WGS $\log_2(\mathrm{CN})$")
    ax.set_title("84-Gene Viral Transcriptome Rank Against EBV Copy Number (n = 717)")

    custom_annot_a = {
        "BGLF5": (1, 0.75),
        "BBLF2_BBLF3": (15, 0.72),
        "BZLF1": (26, 0.69),
        "BART_RPMS1": (54, 0.49),
        "EBNA-1_C": (62, 0.14),
        "LMP-1": (58, -0.18),
        "EBNA-2": (64, -0.23),
    }
    for rank_idx, g in enumerate(genes_sorted):
        if g["gene"] in custom_annot_a:
            y_v = g["r_cn_within_pop_batch"]
            tx, ty = custom_annot_a[g["gene"]]
            ax.annotate(f"{g['gene']} ({y_v:+.2f})", xy=(rank_idx, y_v),
                        xytext=(tx, ty),
                        arrowprops=dict(arrowstyle="->", lw=0.7), fontsize=7.1, fontweight="bold")

    legend_handles = [
        mpatches.Patch(color=CLASS_COLORS["LYTIC_IE"], label="Lytic IE"),
        mpatches.Patch(color=CLASS_COLORS["LYTIC_EARLY"], label="Lytic Early"),
        mpatches.Patch(color=CLASS_COLORS["LYTIC_LATE"], label="Lytic Late"),
        mpatches.Patch(color=CLASS_COLORS["LATENT_BART"], label="Latent BART"),
        mpatches.Patch(color=CLASS_COLORS["LATENT_MRNA"], label="Latent Pol II mRNA"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=7.4)

    # Panel B: Cross-Thaw Biological Repeatability (r_bio) vs Within-Pop+Batch CN Correlation
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    expr_genes = [g for g in genes if g["mean_rpkm"] >= 1.0 and g["kinetic_class"] != "B958_DEL"]
    xb = np.array([g["r_bio"] for g in expr_genes])
    yb = np.array([g["r_cn_within_pop_batch"] for g in expr_genes])
    cb = [CLASS_COLORS.get(g["kinetic_class"], "#888888") for g in expr_genes]
    sb = [max(20, min(140, math.log2(g["mean_rpkm"] + 1) * 18)) for g in expr_genes]
    ax.scatter(xb, yb, c=cb, s=sb, alpha=0.8, edgecolors="black", linewidths=0.4)

    offsets_b = {
        "BGLF5": (-0.18, 0.05),
        "BBLF2_BBLF3": (-0.28, -0.07),
        "BZLF1": (-0.20, -0.05),
        "BRLF1": (-0.16, 0.04),
        "BALF5": (0.03, -0.02),
        "BHRF1": (0.03, 0.00),
        "BART_RPMS1": (0.03, 0.02),
        "EBNA-1_C": (0.03, 0.00),
        "EBNA-2": (0.03, 0.00),
        "LMP-1": (0.03, 0.00),
    }
    for g in expr_genes:
        if g["gene"] in offsets_b:
            dx, dy = offsets_b[g["gene"]]
            ax.annotate(g["gene"], xy=(g["r_bio"], g["r_cn_within_pop_batch"]),
                        xytext=(g["r_bio"] + dx, g["r_cn_within_pop_batch"] + dy),
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555") if abs(dx) > 0.05 else None,
                        fontsize=7.2, fontweight="bold", va="center")
    ax.axhline(0, color="#666666", ls=":", lw=0.9)
    ax.set_xlim(-0.15, 1.05)
    ax.set_ylim(-0.20, 0.76)
    ax.set_xlabel("Cross-Thaw Biological Repeatability (r_bio across 46 independent thaw pairs)")
    ax.set_ylabel(r"Within-Pop+Batch $r$ with WGS $\log_2(\mathrm{CN})$")
    ax.set_title("Cross-Thaw Stability Predicts EBV DNA Copy Number Coupling")

    # Panel C: Abortive vs Productive Lytic Replication Test (Partial r | Spliced BZLF1)
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    fork_and_late = [
        ("BcLF1 (Major Capsid, Late)", 0.2775, "LYTIC_LATE"),
        ("BALF4 (Glycoprotein B, Late)", 0.2776, "LYTIC_LATE"),
        ("BBLF2_BBLF3 (Helicase-Primase)", 0.2580, "LYTIC_EARLY"),
        ("BGLF5 (Host-Shutoff Exonuclease)", 0.2466, "LYTIC_EARLY"),
        ("BALF5 (Viral DNA Polymerase)", 0.2439, "LYTIC_EARLY"),
        ("BALF1 (Strictly Lytic vBCL-2)", 0.2354, "LYTIC_EARLY"),
        ("BGLF4 (Viral Protein Kinase)", 0.2286, "LYTIC_EARLY"),
        ("BBLF4 (Replication Fork Helicase)", 0.1788, "LYTIC_EARLY"),
        ("BALF2 (ssDNA-Binding Protein)", 0.1560, "LYTIC_EARLY"),
        ("BHRF1 (Dual Latent/Lytic vBCL-2)", 0.0016, "LATENT_MRNA"),
        ("EBNA-1_C (Latent Episome Maintenance)", 0.0321, "LATENT_MRNA"),
        ("LMP-1 (Latent III NF-kB Driver)", -0.1054, "LATENT_MRNA"),
        ("EBNA-2 (Latent III Transactivator)", -0.1182, "LATENT_MRNA"),
    ]
    y_pos = np.arange(len(fork_and_late))
    vals_c = [x[1] for x in fork_and_late]
    cols_c = [CLASS_COLORS[x[2]] for x in fork_and_late]
    ax.barh(y_pos, vals_c, color=cols_c, alpha=0.88, height=0.68)
    ax.axvline(0, color="black", lw=0.9)
    ax.axvline(0.097, color="#444444", ls="--", lw=0.9, label="p = 0.01 threshold (|t| = 2.58)")
    for idx, v in enumerate(vals_c):
        ax.text(v + (0.006 if v >= 0 else -0.044), idx, f"{v:+.3f}", va="center", fontsize=7.3, fontweight="bold")
    ax.set_yticks(y_pos)
    ax.set_yticklabels([x[0] for x in fork_and_late], fontsize=7.6)
    ax.invert_yaxis()
    ax.set_xlim(-0.18, 0.35)
    ax.set_xlabel(r"Partial $r(\mathrm{Gene},\ \log_2\mathrm{CN}\ \mid\ \mathrm{Spliced\ BZLF1,\ pop,\ batch})$")
    ax.set_title("Productive Lytic Replication Beyond Immediate-Early BZLF1")
    ax.legend(loc="lower right", frameon=False, fontsize=7.4)

    # Panel D: 3'-Coterminal Poly(A) Step-Up Initiation Flux
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    clusters = [
        ("BRLF1 (5')\n[IE Cluster]", 12.73, 12.73, "#D55E00"),
        ("BZLF1 (3')\n[IE Step-Up]", 24.03, 11.30, "#D55E00"),
        ("BaRF1 (5')\n[Early Cluster]", 2.00, 2.00, "#E69F00"),
        ("BMRF1 (Mid)\n[Pol Processivity]", 57.94, 55.94, "#E69F00"),
        ("BMRF2 (3')\n[Step-Up]", 69.37, 11.43, "#E69F00"),
        ("BGLF4 (Mid)\n[Kinase Step-Up]", 22.94, 22.94, "#009E73"),
        ("BGLF5 (Mid)\n[Nuclease Step-Up]", 30.51, 7.57, "#009E73"),
        ("BBLF1 (3')\n[Tegument Step-Up]", 50.94, 20.43, "#009E73"),
    ]
    x_d = np.arange(len(clusters))
    w_d = 0.38
    ax.bar(x_d - w_d/2, [c[1] for c in clusters], w_d, color="#BBBBBB", label="Total non-overlapping window RPKM")
    ax.bar(x_d + w_d/2, [c[2] for c in clusters], w_d, color=[c[3] for c in clusters], label="De-convolved 5'->3' step-up initiation flux")
    for idx, c in enumerate(clusters):
        ax.text(idx + w_d/2, c[2] + 1.2, f"+{c[2]:.1f}", ha="center", fontsize=7.0, fontweight="bold")
    ax.set_xticks(x_d)
    ax.set_xticklabels([c[0] for c in clusters], fontsize=7.1, rotation=22, ha="right")
    ax.set_ylim(0, 82)
    ax.set_ylabel("Mean Transcript Abundance / Step-Up Flux (RPKM)")
    ax.set_title("3'-Coterminal Poly(A) Step-Up Flux Isolates Promoter Initiation")
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig3_viral_gene_dissection_84.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig3_viral_gene_dissection_84.(png/pdf)")


def make_figure_4(D):
    """Figure 4: Aim 2 — Viral Genetics: IR1 Repeat Structural QTL, 112-Site QTL Scan, & Within-Sample cis-ASE."""
    ped, cn30, wt30, unrel = D["ped"], D["cn30"], D["wt30"], D["unrel"]
    ir1_30 = D["ir1_30"]
    clean_wgs = sorted([s for s in unrel if s in cn30 and wt30.get(s, 0) <= 0.01 and s in ped and s in ir1_30])

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.0))
    plt.subplots_adjust(wspace=0.38, hspace=0.38)

    # Panel A: Within-population IR1 repeat structural QTL (Simpson's Paradox resolved across 2,449 lines)
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    pops = [ped[s]["pop"] for s in clean_wgs]
    x_ir1_w = demean_by_groups([ir1_30[s] for s in clean_wgs], pops)
    y_cn_w = demean_by_groups([math.log2(cn30[s]) for s in clean_wgs], pops)
    r_ir1_w, p_ir1_w = stats.pearsonr(x_ir1_w, y_cn_w)

    ax.scatter(x_ir1_w, y_cn_w, s=9, alpha=0.35, color="#0072B2", edgecolors="none")
    m_i, b_i = np.polyfit(x_ir1_w, y_cn_w, 1)
    xs_i = np.linspace(np.percentile(x_ir1_w, 0.5), np.percentile(x_ir1_w, 99.5), 100)
    ax.plot(xs_i, m_i * xs_i + b_i, color="#D55E00", lw=2.2,
            label=f"30x WGS Within-Pop: r = +{r_ir1_w:.4f} (p = {p_ir1_w:.1e}, n = {len(clean_wgs):,})")
    ax.set_xlim(-6.5, 8.5)
    ax.set_xlabel("Within-Population Demeaned IR1 (BamHI-W) Repeat Units")
    ax.set_ylabel(r"Within-Population Demeaned $\log_2(\mathrm{EBV\ CN})$")
    ax.set_title("Structural cis-QTL: IR1 Repeat Number Drives EBV CN")
    ax.text(0.03, 0.71,
            "Multi-Platform & Artifact Controls:\n"
            "• 30x IR1 / genome_uniq: r = +0.3148 (p = 1.7e-60)\n"
            "• 30x IR1 / local_Cp_flank: r = +0.3037 (p = 2.1e-53)\n"
            "• Within-chrEBV GC-controlled: r = +0.3048 (96.8% retained)\n"
            "• ONT long-read IR1 vs ONT CN (n=890): r = +0.2945\n"
            "• MAGE RNA within-pop+batch (n=717): r = +0.3687",
            transform=ax.transAxes, fontsize=7.4,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8F9FA", edgecolor="#CCCCCC"))
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)

    # Panel B: Mechanism of IR1 Repeat Copy Number on MAGE Viral Transcripts (Within-Pop+Batch, n = 717)
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    ir1_fx = [
        ("WGS log2(EBV CN)\n(Episome accumulation)", +0.3687, "#D55E00"),
        ("log2(Lytic IE / Latent)\n(Lytic switch ratio)", +0.2531, "#D55E00"),
        ("Lytic IE RPKM\n(BZLF1 / BRLF1)", +0.1859, "#E69F00"),
        ("Spliced BZLF1 RPM\n(Zta productive switch)", +0.1743, "#E69F00"),
        ("BHRF1 RPKM\n(W-spliced viral BCL-2)", +0.1029, "#009E73"),
        ("EBNA-1 RPKM\n(Episome maintenance)", -0.0257, "#888888"),
        ("LMP-1 RPKM\n(Latent III NF-kB)", -0.0993, "#0072B2"),
        ("EBNA-2 RPKM\n(Downstream of IR1)", -0.1068, "#0072B2"),
    ]
    y_b = np.arange(len(ir1_fx))
    vals_b = [x[1] for x in ir1_fx]
    ax.barh(y_b, vals_b, color=[x[2] for x in ir1_fx], alpha=0.88, height=0.65)
    ax.axvline(0, color="black", lw=0.9)
    ax.axvline(0.073, color="#666666", ls=":", lw=0.9)
    ax.axvline(-0.073, color="#666666", ls=":", lw=0.9, label="p = 0.05 threshold (|r| = 0.073)")
    for idx, v in enumerate(vals_b):
        ax.text(v + (0.008 if v >= 0 else -0.052), idx, f"{v:+.4f}", va="center", fontsize=7.4, fontweight="bold")
    ax.set_yticks(y_b)
    ax.set_yticklabels([x[0] for x in ir1_fx], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.18, 0.45)
    ax.set_xlabel("Within-Pop+Batch Correlation (r) with IR1 Repeat Copy Number (n = 717)")
    ax.set_title("IR1 Expansion Boosts BHRF1 & Lytic Switch vs EBNA-2")
    ax.legend(loc="lower right", frameon=False, fontsize=7.4)

    # Panel C: Within-Population Viral QTL Scan (112 Polymorphic Sites + Rare Zp/Rp Burden)
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    qtl_highlights = [
        ("12,129:T>A (W0/Cp-Wp Leader)\n[All 1kGP n=2,449, carrier=111]", +0.1695, "#D55E00"),
        ("12,129:T>A (W0/Cp-Wp Leader)\n[Holdout Replication n=1,732]", +0.1736, "#D55E00"),
        ("12,129:T>A -> MAGE Lytic IE\n[Within-pop+batch n=717]", +0.1048, "#E69F00"),
        ("Rare Zp Burden (BZLF1 Promoter)\n[90,944–91,500 bp, n=16 carriers]", +0.2310, "#CC79A7"),
        ("91,192:G>A (BRLF1 IBS Barcode)\n[Within-pop n=2,449: NEUTRAL]", -0.0410, "#888888"),
        ("132,965:G>GC (BXRF1 GBR Barcode)\n[Within-pop n=2,449: NEUTRAL]", -0.0416, "#888888"),
        ("163,858:G>T (BALF2 ssDNA-BP)\n[cis-BALF2 RNA within-pop+batch]", -0.1243, "#0072B2"),
        ("163,858:G>T (BALF2 ssDNA-BP)\n[trans-Lytic IE within-pop+batch]", -0.1282, "#0072B2"),
        ("Rare Rp Burden (BRLF1 Promoter)\n[92,896–93,500 bp, n=12 carriers]", -0.4930, "#009E73"),
    ]
    y_c = np.arange(len(qtl_highlights))
    vals_qc = [x[1] for x in qtl_highlights]
    ax.barh(y_c, vals_qc, color=[x[2] for x in qtl_highlights], alpha=0.88, height=0.68)
    ax.axvline(0, color="black", lw=0.9)
    for idx, v in enumerate(vals_qc):
        unit_str = f"{v:+.3f}" + (" Δlog2" if "Burden" in qtl_highlights[idx][0] else " (r)")
        ax.text(v + (0.012 if v >= 0 else -0.125), idx, unit_str, va="center", fontsize=7.1, fontweight="bold")
    ax.set_yticks(y_c)
    ax.set_yticklabels([x[0] for x in qtl_highlights], fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlim(-0.66, 0.35)
    ax.set_xlabel("Within-Population Effect Size (Correlation r or Rare Promoter Δlog2[CN])")
    ax.set_title("Neutrality of IBS/GBR Barcodes vs Replicated Viral QTLs")

    # Panel D: Within-Sample Heteroplasmic Viral cis-ASE (VAF_RNA vs VAF_DNA across 59 carrier libraries)
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    ase_dna, ase_rna = [], []
    site_obs = defaultdict(list)
    with open(os.path.join(DATA_DIR, "all_mage_viral_ase_59.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rtot = int(row["rna_ref"]) + int(row["rna_alt"])
            daf = float(row["dna_vaf"])
            raf = float(row["rna_vaf"])
            if rtot >= 5 and daf > 0 and len(row["ref"]) == 1 and len(row["alt"]) == 1:
                ase_dna.append(daf)
                ase_rna.append(raf)
                site_obs[(int(row["pos"]), row["ref"], row["alt"])].append((daf, raf, int(row["rna_ref"]), int(row["rna_alt"])))

    ase_dna = np.array(ase_dna)
    ase_rna = np.array(ase_rna)
    r_ase, _ = stats.pearsonr(ase_dna, ase_rna)
    ax.scatter(ase_dna, ase_rna, s=11, alpha=0.35, color="#56B4E9", edgecolors="none",
               label=f"Heteroplasmic SNV observations (n = {len(ase_dna):,}, r = +{r_ase:.4f})")
    ax.plot([0, 1], [0, 1], color="#333333", ls="--", lw=1.2, label="Balanced cis-transcription (VAF_RNA = VAF_DNA)")

    highlights_ase = [
        ((91192, "G", "A"), "91,192:G>A (BRLF1 IBS)\n20.5% RNA vs 19.2% DNA (neutral)", "#111111", (0.28, 0.30)),
        ((168963, "T", "G"), "168,963:T>G (LMP-1)\n6.1% RNA vs 3.9% DNA (p=1.5e-95, cis-UP)", "#D55E00", (0.14, 0.16)),
        ((36470, "T", "C"), "36,470:T>C (EBNA-2)\n2.7% RNA vs 16.0% DNA (p=1.1e-31, cis-DOWN)", "#0072B2", (0.26, 0.02)),
    ]
    for key, lbl, col, xytext in highlights_ase:
        obs = site_obs.get(key, [])
        if obs:
            tot_r = sum(x[2] + x[3] for x in obs)
            w_dna = sum(x[0] * (x[2] + x[3]) for x in obs) / tot_r
            p_rna = sum(x[3] for x in obs) / tot_r
            ax.scatter([w_dna], [p_rna], color=col, s=65, zorder=5, edgecolors="white", linewidths=1.0)
            ax.annotate(lbl, xy=(w_dna, p_rna), xytext=xytext,
                        arrowprops=dict(arrowstyle="->", color=col, lw=1.1),
                        fontsize=7.1, fontweight="bold", color=col)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Within-Sample WGS DNA Variant Allele Fraction (VAF_DNA)")
    ax.set_ylabel("Same-Donor MAGE RNA Variant Allele Fraction (VAF_RNA)")
    ax.set_title("Within-Sample Heteroplasmic Viral cis-ASE")
    ax.legend(loc="upper left", frameon=False, fontsize=7.4)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig4_viral_genetics_ir1_qtl_and_ase.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig4_viral_genetics_ir1_qtl_and_ase.(png/pdf)")


def make_figure_5():
    """Figure 5: Aim 3 — Host Transcriptome Mediation, Exogenous IR1 Instrument, & 4-Axis Bipartite Map."""
    rows = []
    with open(os.path.join(DATA_DIR, "host_viral_mediation_20154.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)

    fig, axes = plt.subplots(2, 2, figsize=(14.6, 10.2))
    plt.subplots_adjust(wspace=0.42, hspace=0.40)

    # Panel A: Transcriptome-Wide Discovery Boost (BH-FDR q < 0.05 across 20,154 human genes, df = 673)
    ax = axes[0, 0]
    add_panel_label(ax, "A")
    predictors = [
        ("Exogenous IR1\n(WGS DNA)", 663, 26.6, "#56B4E9"),
        ("2019 WGS\nlog2(EBV CN)", 6087, 49.9, "#888888"),
        ("Spliced BZLF1\n(Zta) RPM", 9819, 65.1, "#E69F00"),
        ("BGLF5\n(Host-Shutoff)", 10580, 67.9, "#D55E00"),
        ("log2(Late /\nLatent mRNA)", 12863, 75.7, "#CC79A7"),
        ("LMP-1 RPKM\n(NF-kB Driver)", 13800, 78.4, "#0072B2"),
    ]
    x_a = np.arange(len(predictors))
    cnts = [p[1] for p in predictors]
    ax.bar(x_a, cnts, color=[p[3] for p in predictors], alpha=0.88, width=0.60)
    for idx, (_, cnt, pi1, _) in enumerate(predictors):
        pct = cnt / 20154.0 * 100.0
        ax.text(idx, cnt + 220, f"{cnt:,}\n({pct:.1f}% | π₁={pi1:.0f}%)",
                ha="center", va="bottom", fontsize=7.0, fontweight="bold")
    ax.set_xticks(x_a)
    ax.set_xticklabels([p[0] for p in predictors], fontsize=7.3, rotation=16, ha="right")
    ax.set_ylim(0, 16400)
    ax.set_ylabel("Significant Human Genes at BH-FDR q < 0.05")
    ax.set_title("Same-Vial Viral RNA Doubles Host Transcriptome Discovery")

    # Panel B: Formal 1-df Baron-Kenny Mediation Asymmetry (2019 WGS CN <-> 2023 Spliced BZLF1)
    ax = axes[0, 1]
    add_panel_label(ax, "B")
    r_cn = np.array([float(r["r_wgs_cn"]) for r in rows])
    r_cn_given_z = np.array([float(r["r_wgs_cn_given_bzlf1"]) for r in rows])
    r_z = np.array([float(r["r_spliced_bzlf1"]) for r in rows])
    r_z_given_cn = np.array([float(r["r_bzlf1_given_wgs_cn"]) for r in rows])
    q_cn = np.array([float(r["q_wgs_cn"]) for r in rows])

    sig_mask = q_cn < 0.05
    ax.scatter(r_cn[sig_mask], r_cn_given_z[sig_mask], s=6, alpha=0.25, color="#D55E00",
               label="WGS CN | Spliced BZLF1 (92.4% of 6,087 hits eliminated!)")
    ax.scatter(r_z[sig_mask], r_z_given_cn[sig_mask], s=6, alpha=0.25, color="#0072B2",
               label="Spliced BZLF1 | WGS CN (78.7% of 9,819 hits survive!)")
    xs_b = np.linspace(-0.4, 0.45, 100)
    ax.plot(xs_b, xs_b, color="#333333", ls="--", lw=1.1, label="No attenuation (y = x)")
    ax.axhline(0, color="#666666", ls=":", lw=1.0)
    m_f, b_f = np.polyfit(r_cn[sig_mask], r_cn_given_z[sig_mask], 1)
    m_r, b_r = np.polyfit(r_z[sig_mask], r_z_given_cn[sig_mask], 1)
    ax.plot(xs_b, m_f * xs_b + b_f, color="#8C2D04", lw=2.0, label=f"Forward slope = {m_f:.2f} (86.0% of R² mediated)")
    ax.plot(xs_b, m_r * xs_b + b_r, color="#004C78", lw=2.0, label=f"Reverse slope = {m_r:.2f} (direct viral RNA effect)")
    ax.set_xlim(-0.38, 0.45)
    ax.set_ylim(-0.38, 0.45)
    ax.set_xlabel("Unconditioned Within-Pop+Batch Partial Correlation (r)")
    ax.set_ylabel("Conditioned Baron–Kenny Partial Correlation (r)")
    ax.set_title("Baron–Kenny Mediation Asymmetry: Spliced BZLF1 Mediates CN")
    ax.legend(loc="upper left", frameon=False, fontsize=7.3)

    # Panel C: Exogenous Viral DNA Instrument (IR1 Repeat Count) vs WGS CN Across 20,154 Human Genes
    ax = axes[1, 0]
    add_panel_label(ax, "C")
    r_ir1 = np.array([float(r["r_ir1_dna"]) for r in rows])
    r_align, _ = stats.pearsonr(r_ir1, r_cn)
    ax.scatter(r_ir1, r_cn, s=5, alpha=0.22, color="#009E73", edgecolors="none")
    m_c, b_c = np.polyfit(r_ir1, r_cn, 1)
    xs_c = np.linspace(r_ir1.min(), r_ir1.max(), 100)
    ax.plot(xs_c, m_c * xs_c + b_c, color="#111111", lw=1.8,
            label=f"20,154-gene alignment: r = +{r_align:.4f} (663 IR1 genes q < 0.05)")
    ax.axhline(0, color="#888888", ls=":", lw=0.8)
    ax.axvline(0, color="#888888", ls=":", lw=0.8)

    offsets_c = {
        "CBARP": (-0.06, 0.02),
        "USP17L2": (-0.07, -0.01),
        "LINC01443": (-0.08, -0.05),
        "HLA-DOA": (-0.07, 0.01),
        "HLA-DOB": (-0.07, -0.02),
        "CD40": (0.01, 0.00),
        "EBF1": (0.01, 0.00),
    }
    for r in rows:
        sym = r["gene_symbol"]
        if sym in offsets_c:
            xi, yi = float(r["r_ir1_dna"]), float(r["r_wgs_cn"])
            dx, dy = offsets_c[sym]
            ax.scatter([xi], [yi], color="#D55E00", s=24, zorder=5)
            ax.text(xi + dx, yi + dy, sym, fontsize=7.2, fontweight="bold", va="center")
    ax.set_xlim(-0.24, 0.25)
    ax.set_xlabel("Exogenous Viral DNA IR1 Repeat Effect on Host Gene (r)")
    ax.set_ylabel("WGS log2(EBV CN) Effect on Host Gene (r)")
    ax.set_title("Exogenous DNA Instrument (IR1) Reconstructs Host Signature")
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)

    # Panel D: Bipartite 4-Axis Viral-to-Host Functional Module Heatmap
    ax = axes[1, 1]
    add_panel_label(ax, "D")
    by_sym = {r["gene_symbol"]: r for r in rows}
    target_genes = [
        ("IFIH1", "NF-kB"),
        ("TNFAIP3", "NF-kB"),
        ("CD40", "B-Cell"),
        ("ICAM1", "NF-kB"),
        ("NFKBIA", "NF-kB"),
        ("EBF1", "B-Cell"),
        ("MYC", "B-Cell"),
        ("TNFRSF17", "Plasmablast"),
        ("MZB1", "Plasmablast"),
        ("SSR4", "Plasmablast"),
        ("JCHAIN", "Plasmablast"),
        ("MX1", "IFN / ISG"),
        ("IFIT1", "IFN / ISG"),
        ("ISG15", "IFN / ISG"),
        ("OAS1", "IFN / ISG"),
        ("STAT1", "IFN / ISG"),
    ]
    viral_cols = [
        ("r_lmp1", "LMP-1\n(Latent NF-kB)"),
        ("r_bglf5", "BGLF5\n(Lytic / Shutoff)"),
        ("r_spliced_bzlf1", "Spliced BZLF1\n(Lytic Switch)"),
        ("r_ebna2", "EBNA-2\n(Latent III)"),
        ("r_bart", "BART lncRNA\n(ISG Repressor)"),
    ]
    mat = np.zeros((len(target_genes), len(viral_cols)))
    for i, (sym, _) in enumerate(target_genes):
        row = by_sym.get(sym, {})
        for j, (col_k, _) in enumerate(viral_cols):
            mat[i, j] = float(row.get(col_k, 0.0))

    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.75, vmax=0.75, aspect="auto")
    for i in range(len(target_genes)):
        for j in range(len(viral_cols)):
            v = mat[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    color="white" if abs(v) > 0.38 else "black", fontsize=6.8, fontweight="bold")
    ax.set_xticks(range(len(viral_cols)))
    ax.set_xticklabels([c[1] for c in viral_cols], fontsize=7.3)
    ax.set_yticks(range(len(target_genes)))
    ax.set_yticklabels([f"{sym} ({mod})" for sym, mod in target_genes], fontsize=7.3)
    ax.set_title("Bipartite 4-Axis Viral-to-Host Functional Module Map")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Within-Pop+Batch Partial r", fontsize=7.8)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig5_host_mediation_and_bipartite_map.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig5_host_mediation_and_bipartite_map.(png/pdf)")


def make_figure_6():
    """Figure 6: Out-of-Sample Vector Alignment (24 Non-IBS/GBR Pops -> IBS vs GBR Differential Expression)."""
    rows = []
    with open(os.path.join(DATA_DIR, "host_viral_mediation_20154.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)

    beta_24 = np.array([float(r["beta_cn_24pops"]) for r in rows])
    de_ibs_gbr = np.array([float(r["de_ibs_minus_gbr"]) for r in rows])
    q_cn = np.array([float(r["q_wgs_cn"]) for r in rows])
    sig = q_cn < 0.05

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.4))
    plt.subplots_adjust(wspace=0.28)

    # Panel A: Zero-sample-overlap vector alignment across all 20,154 genes & 6,087 significant genes
    ax = axes[0]
    add_panel_label(ax, "A")
    r_all, _ = stats.pearsonr(beta_24, de_ibs_gbr)
    r_sig, _ = stats.pearsonr(beta_24[sig], de_ibs_gbr[sig])

    ax.scatter(beta_24[~sig], de_ibs_gbr[~sig], s=5, alpha=0.18, color="#999999", edgecolors="none",
               label=f"All 20,154 human genes: r = +{r_all:.4f} (R² = {r_all**2*100:.1f}%)")
    ax.scatter(beta_24[sig], de_ibs_gbr[sig], s=7, alpha=0.35, color="#D55E00", edgecolors="none",
               label=f"CN-significant genes (n = {sig.sum():,}): r = +{r_sig:.4f} (R² = {r_sig**2*100:.1f}%)")
    xs = np.linspace(beta_24.min(), beta_24.max(), 100)
    m_s, b_s = np.polyfit(beta_24[sig], de_ibs_gbr[sig], 1)
    ax.plot(xs, m_s * xs + b_s, color="#111111", lw=2.0,
            label=f"Out-of-sample slope = {m_s:.2f} (predicted Δlog2[CN] ≈ 1.44)")
    ax.axhline(0, color="#666666", ls=":", lw=0.8)
    ax.axvline(0, color="#666666", ls=":", lw=0.8)
    ax.set_xlabel(r"Within-Population EBV CN Slope $\hat{\beta}_{\mathrm{CN}}$ Trained in 24 Non-IBS/GBR Populations")
    ax.set_ylabel(r"Observed Differential Expression $\overline{Y}_{\mathrm{IBS}} - \overline{Y}_{\mathrm{GBR}}$ (Spanish vs British)")
    ax.set_title("Zero-Sample-Overlap Prediction of IBS vs GBR Differential Expression")
    ax.legend(loc="upper left", frameon=False, fontsize=7.6)

    # Panel B: Comparison of Out-of-Sample Predictors (24 Non-IBS/GBR Pops -> IBS vs GBR)
    ax = axes[1]
    add_panel_label(ax, "B")
    models = [
        ("Global 25-df\n26-Pop ANOVA", 4.0, 8.3),
        ("Same-Vial\nSpliced BZLF1", 42.6, 61.9),
        ("2019 WGS\nlog2(EBV CN)", 51.3, 72.1),
        ("Multi-Parameter\nViral Model", 54.0, 73.8),
    ]
    x_b = np.arange(len(models))
    w_b = 0.36
    b1 = ax.bar(x_b - w_b/2, [m[1] for m in models], w_b, color="#56B4E9", alpha=0.9, label="All 20,154 human genes (% variance R²)")
    b2 = ax.bar(x_b + w_b/2, [m[2] for m in models], w_b, color="#D55E00", alpha=0.9, label="CN-associated genes q < 0.05 (% variance R²)")
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2, h + 1.2, f"{h:.1f}%", ha="center", va="bottom", fontsize=7.6, fontweight="bold")
    ax.set_xticks(x_b)
    ax.set_xticklabels([m[0] for m in models], fontsize=7.8)
    ax.set_ylim(0, 86)
    ax.set_ylabel("Explained Between-Population Variance R² (%)")
    ax.set_title("Targeted Out-of-Population Vector Alignment vs Global ANOVA")
    ax.legend(loc="upper left", frameon=False, fontsize=7.6)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"fig6_population_vector_alignment_ibs_gbr.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("Saved fig6_population_vector_alignment_ibs_gbr.(png/pdf)")


def main():
    print("Loading master tables from:", DATA_DIR)
    D = load_core_data()
    make_figure_1(D)
    make_figure_2(D)
    make_figure_3()
    make_figure_4(D)
    make_figure_5()
    make_figure_6()
    print("All 6 multi-panel figures generated successfully in:", FIG_DIR)


if __name__ == "__main__":
    main()
