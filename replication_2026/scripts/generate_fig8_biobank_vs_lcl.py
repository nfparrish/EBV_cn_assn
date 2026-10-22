#!/usr/bin/env python3
"""
Generates Figure 8: Cross-Study Architecture of In Vivo Biobank EBV GWAS
(UK Biobank + All of Us Blood & Saliva; Lareau/Kamitaki/Okada 2026) vs.
Viral-Conditioned B-Cell-Intrinsic LCL Host & Viral Genetics (1kGP + MAGE).
"""
import csv
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
FIG_DIR = os.path.join(ROOT, "figures")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(16.5, 12.2), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.05, 1.0], wspace=0.26, hspace=0.34)

    # =========================================================================
    # PANEL A: High-Prior Biobank vs. LCL Locus Comparison Across Compartments
    # =========================================================================
    axA = fig.add_subplot(gs[0, 0])
    loci_rows = []
    with open(os.path.join(DATA_DIR, "biobank_ebv_loci_comparison.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            loci_rows.append(r)

    # Select 22 representative headline loci across the 4 compartments
    selected_genes = [
        # MHC Class I/II (4)
        "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1",
        # ER Trimming & T/NK Cell-Extrinsic Immune (8)
        "ERAP2", "ERAP1", "PTPN22", "LY9", "SH2B3", "CTLA4", "GSDMB", "CD70", "CD27",
        # Shared B-Cell Intrinsic / Entry Co-Receptor TFs (5)
        "IRF4", "EBF1", "CIITA", "CD40", "PAX5",
        # Epigenetic CpG Methylation & Heterochromatin Silencers (4)
        "DNMT1", "TET2", "ZNF814", "DNMT3B",
    ]
    row_by_g = {r["gene"]: r for r in loci_rows}
    plot_loci = [row_by_g[g] for g in selected_genes if g in row_by_g]

    def safe_neglog10(p_str, cap=25.0):
        try:
            pv = float(p_str)
            if pv <= 0.0:
                return cap
            return min(cap, -math.log10(pv))
        except Exception:
            return 0.0

    y_pos = np.arange(len(plot_loci))
    blood_nlp = [safe_neglog10(r["biobank_blood_p"], cap=25.0) for r in plot_loci]
    saliva_nlp = [safe_neglog10(r["biobank_saliva_p"], cap=25.0) for r in plot_loci]
    # Exact Biobank lead SNP in LCL WGS_CN_ViralCond
    lcl_exact_nlp = [safe_neglog10(r["lcl_exact_blood_lead_CN_p"], cap=25.0) for r in plot_loci]
    # LCL Window min P across ViralCond phenotypes (min of WGS_CN, RNA_IELatRatio, RNA_ProdFork, RNA_LMP1)
    lcl_win_nlp = [
        max(
            safe_neglog10(r["lcl_window_minP_WGS_CN_ViralCond"], cap=25.0),
            safe_neglog10(r["lcl_window_minP_RNA_IELatRatio_ViralCond"], cap=25.0),
            safe_neglog10(r["lcl_window_minP_RNA_ProdFork_ViralCond"], cap=25.0),
            safe_neglog10(r["lcl_window_minP_RNA_LMP1_ViralCond"], cap=25.0),
        )
        for r in plot_loci
    ]

    bar_h = 0.20
    axA.barh(y_pos - 1.5 * bar_h, blood_nlp, height=bar_h, color="#1f4e78", alpha=0.92, label="Biobank Blood ($N\\approx 900\\mathrm{k}$, capped at 25)")
    axA.barh(y_pos - 0.5 * bar_h, saliva_nlp, height=bar_h, color="#2e8b57", alpha=0.88, label="Biobank Saliva Shedding ($N\\approx 245\\mathrm{k}$)")
    axA.barh(y_pos + 0.5 * bar_h, lcl_exact_nlp, height=bar_h, color="#999999", alpha=0.90, label="LCL Exact Biobank Lead SNP (`WGS_CN_ViralCond`)")
    axA.barh(y_pos + 1.5 * bar_h, lcl_win_nlp, height=bar_h, color="#c0392b", alpha=0.92, label="LCL Window Lead ($\\pm 100$ kb, Viral-Conditioned)")

    axA.axvline(-math.log10(5e-8), color="#1f4e78", linestyle="--", linewidth=1.0, alpha=0.7)
    axA.axvline(-math.log10(0.05), color="#555555", linestyle=":", linewidth=0.9, alpha=0.7)
    axA.text(-math.log10(5e-8) + 0.3, len(plot_loci) - 0.6, "$P=5\\times 10^{-8}$", fontsize=7.5, color="#1f4e78")

    # Compartment dividers
    for sep_y in (3.5, 12.5, 17.5):
        axA.axhline(sep_y, color="#cccccc", linestyle="-", linewidth=0.9)

    txt_box = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.88)
    axA.text(14.2, 2.0, "MHC Class I/II ($P_{\\mathrm{Blood}}<10^{-300}$): Exact lead SNPs null in LCLs", fontsize=7.6, color="#1f4e78", fontweight="bold", bbox=txt_box)
    axA.text(13.2, 10.6, "Cell-Extrinsic T/NK & ERAP1/2: 100% null at exact\nlead SNPs in T-cell-free B-cell monoculture", fontsize=7.6, color="#2c3e50", fontweight="bold", bbox=txt_box)
    axA.text(6.0, 15.1, "Shared Sub-Threshold\nB-Cell TFs (EBF1, IRF4,\nCIITA, CD40, PAX5)", fontsize=7.4, color="#9c640c", fontweight="bold", bbox=txt_box)
    axA.text(8.2, 19.6, "B-Cell Epigenetic Silencers Unmasked in LCLs:\nDNMT3B ($q=0.022$), ZNF814 ($q=0.085$)", fontsize=7.6, color="#922b21", fontweight="bold", bbox=txt_box)

    axA.set_yticks(y_pos)
    axA.set_yticklabels([r["gene"] for r in plot_loci], fontweight="bold", fontsize=8.3)
    axA.invert_yaxis()
    axA.set_xlabel("Association Strength ($-\\log_{10} P$, capped at 25 for MHC/ERAP)", fontweight="bold")
    axA.set_xlim(0, 28.0)
    axA.set_title("A. High-Prior Biobank Loci vs. Viral-Conditioned LCL Architecture", fontweight="bold", loc="left", pad=8)
    axA.legend(loc="center right", bbox_to_anchor=(0.99, 0.31), frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.95, fontsize=7.0)

    # =========================================================================
    # PANEL B: Unsigned Quantized Block Test (Immune/MHC-Masked Sub-Threshold Tail)
    # =========================================================================
    axB = fig.add_subplot(gs[0, 1])
    uns_rows = []
    with open(os.path.join(DATA_DIR, "biobank_unsigned_block_concordance.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            uns_rows.append(r)

    # Plot LCL block-max chi-square ratio vs background across Biobank tiers for:
    # 1. Blood Unmasked Tier0 vs Masked Tiers 1-4 (WGS_CN_ViralCond & RNA_LMP1_ViralCond)
    # 2. Saliva Masked Tiers 1-4 (RNA_IELatRatio_ViralCond, RNA_BGLF5_ViralCond, WGS_CN_ViralCond)
    tier_Order = [
        ("Tier0_GWSig_or_MHC_P_lt_5e-8", "Tier 0: MHC /\nGWSig Immune\n($P < 5\\times 10^{-8}$)"),
        ("Tier1_SubThreshold_5e-8_to_1e-5", "Tier 1: Strong\nSub-Threshold\n($5\\times 10^{-8}$–$10^{-5}$)"),
        ("Tier2_SubThreshold_1e-5_to_1e-3", "Tier 2: Moderate\nSub-Threshold\n($10^{-5}$–$10^{-3}$)"),
        ("Tier3_Nominal_1e-3_to_0.01", "Tier 3: Nominal\nSub-Threshold\n($10^{-3}$–$10^{-2}$)"),
        ("Tier4_Background_P_ge_0.01", "Tier 4: Null\nBackground\n($P \\geq 0.01$)"),
    ]

    def get_tier_series(cohort, pheno, col="frac_lcl_block_p_lt_1e-3"):
        vals = []
        for t_code, _ in tier_Order:
            mode = "Unmasked_All_Blocks" if t_code.startswith("Tier0") else "Immune_MHC_Masked_SubThreshold"
            match = [float(r[col]) for r in uns_rows if r["biobank_cohort"] == cohort and r["mask_mode"] == mode and r["lcl_phenotype"] == pheno and r["biobank_tier"] == t_code]
            vals.append(match[0] * 100.0 if match else np.nan)
        return vals

    x_t = np.arange(len(tier_Order))
    s_sal_cn = get_tier_series("Saliva_AoU", "WGS_CN_ViralCond")
    s_sal_ie = get_tier_series("Saliva_AoU", "RNA_IELatRatio_ViralCond")
    s_sal_bg = get_tier_series("Saliva_AoU", "RNA_BGLF5_ViralCond")
    s_bld_cn = get_tier_series("Blood_UKB_AoU", "WGS_CN_ViralCond")
    s_bld_pf = get_tier_series("Blood_UKB_AoU", "RNA_ProdFork_ViralCond")

    axB.plot(x_t, s_sal_cn, marker="o", linewidth=2.2, markersize=7, color="#2e8b57", label="Biobank Saliva $\\rightarrow$ LCL `WGS_CN_ViralCond` ($\\rho_{\\mathrm{mean}\\,\\chi^2}=+0.028, P=0.043$)")
    axB.plot(x_t, s_sal_bg, marker="s", linewidth=2.2, markersize=6.5, color="#16a085", linestyle="--", label="Biobank Saliva $\\rightarrow$ LCL `RNA_BGLF5_ViralCond` ($\\rho_{\\mathrm{mean}\\,\\chi^2}=+0.058, P=2.0\\times 10^{-5}$)")
    axB.plot(x_t, s_sal_ie, marker="^", linewidth=2.0, markersize=6.5, color="#27ae60", linestyle=":", label="Biobank Saliva $\\rightarrow$ LCL `RNA_IELatRatio` (Top-10% $\\mathrm{OR}=1.62, P=6.1\\times 10^{-4}$)")
    axB.plot(x_t, s_bld_pf, marker="D", linewidth=2.2, markersize=6.5, color="#c0392b", label="Biobank Blood $\\rightarrow$ LCL `RNA_ProdFork` (Masked Top-10% $\\mathrm{OR}=1.41, P=0.016$)")
    axB.plot(x_t, s_bld_cn, marker="v", linewidth=2.0, markersize=6.5, color="#2980b9", linestyle="-.", label="Biobank Blood $\\rightarrow$ LCL `WGS_CN_ViralCond` (Tier 1 $>$ Tier 0 Immune)")

    # Highlight the Tier 1 Sub-Threshold peak vs Tier 0 Immune peak
    axB.axvspan(0.65, 1.35, color="#f9ebea", alpha=0.55, zorder=0)
    axB.annotate(
        "Sub-Threshold Non-Immune Tier 1\nExceeds Tier 0 MHC/Immune Peaks!\n(e.g., Saliva Tier 1: 30.8% vs 10.1% BG;\nBlood ProdFork Tier 1: 25.8% vs 17.3% Tier 0)",
        xy=(1.0, 26.5),
        xytext=(1.55, 27.2),
        fontsize=8.0,
        fontweight="bold",
        color="#922b21",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff8f6", edgecolor="#c0392b", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.3),
    )

    axB.set_xticks(x_t)
    axB.set_xticklabels([lbl for _, lbl in tier_Order], fontsize=8.0, fontweight="bold")
    axB.set_ylabel("% of 500-kb Genomic Blocks with\nLCL Block-Lead $P < 10^{-3}$", fontweight="bold")
    axB.set_ylim(1.5, 33.8)
    axB.set_title("B. Unsigned Quantized 500-kb Block Test (MHC/Immune Masked)", fontweight="bold", loc="left", pad=8)
    axB.legend(loc="lower left", frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.95, fontsize=7.2)

    # =========================================================================
    # PANEL C: Signed Quantized Block & TWAS Tests (Both Options 2A and 2B)
    # =========================================================================
    axC = fig.add_subplot(gs[1, 0])
    # Show side-by-side:
    # Left group (4 bars): Cross-Study Germline TWAS Correlation r(Z_TWAS_Bio, Z_TWAS_LCL)
    #   1. Unmasked Tier 0 GWSig Immune Peaks (Blood -> WGS_CN_PopAdj: r = -0.254, P = 2.4e-4)
    #   2. Unmasked Tier 0 GWSig Immune Peaks (Blood -> WGS_CN_ViralCond: r = -0.133, P = 0.059)
    #   3. Immune-Masked Sub-Threshold Tier 1 (Blood P < 1e-3 -> WGS_CN_ViralCond: r = +0.535, P = 4.3e-10)
    #   4. Immune-Masked Sub-Threshold Tier 1 (Blood P < 1e-3 -> RNA_BGLF5_ViralCond: r = +0.424, P = 1.7e-6)
    #   5. Immune-Masked Dual-Nominal eQTLs (Blood & LCL P < 0.05 -> RNA_IELatRatio: r = +0.473, P = 0.0031, 64.9% sign)
    #   6. Option 2A Saliva Dual-Nominal 500-kb Blocks (Saliva & LCL P < 0.05 -> RNA_BGLF5: r = +0.143, P = 0.047, 57.0% sign)
    #   7. Option 2B-2 LCL Germline TWAS vs Observed MAGE mRNA (Nominal P<0.05, n=124: r = +0.347, P = 8.0e-5, 66.9% sign)
    #   8. Option 2B-2 LCL Germline TWAS vs Observed MAGE mRNA (PIP>=0.5 & P<0.05, n=62: r = +0.423, P = 6.1e-4, 71.0% sign)

    c_labels = [
        "Unmasked Tier 0 Immune Peaks ($P<5\\times 10^{-8}$)\nBlood $\\leftrightarrow$ LCL `WGS_CN_PopAdj` ($n=204$)",
        "Unmasked Tier 0 Immune Peaks ($P<5\\times 10^{-8}$)\nBlood $\\leftrightarrow$ LCL `WGS_CN_ViralCond` ($n=204$)",
        "Immune-Masked Sub-Threshold Tier 1 ($P<10^{-3}$)\nBlood $\\leftrightarrow$ LCL `WGS_CN_ViralCond` ($n=118$)",
        "Immune-Masked Sub-Threshold Tier 1 ($P<10^{-3}$)\nBlood $\\leftrightarrow$ LCL `RNA_BGLF5_ViralCond` ($n=118$)",
        "Immune-Masked Dual-Nominal eQTLs ($P<0.05$)\nBlood $\\leftrightarrow$ LCL `RNA_IELatRatio` ($n=37$)",
        "Option 2A Dual-Nominal 500-kb Blocks ($P<0.05$)\nSaliva $\\leftrightarrow$ LCL `RNA_BGLF5` ($n=193$ blocks)",
        "Option 2B-2: LCL TWAS vs. MAGE mRNA Coupling\nNominal LCL eQTLs ($P<0.05, n=124$ genes)",
        "Option 2B-2: LCL TWAS vs. MAGE mRNA Coupling\nFine-Mapped $\\mathrm{PIP}\\geq 0.5$ ($P<0.05, n=62$ genes)",
    ]
    c_rvals = [-0.2543, -0.1325, +0.5348, +0.4241, +0.4729, +0.1431, +0.3469, +0.4230]
    c_sign_txt = [
        "41.7% sign ($P=0.021$, $r=-0.25, P=2.4\\times 10^{-4}$)",
        "44.1% sign ($r=-0.13, P=0.059$)",
        "54.2% sign ($r=+0.53, P=4.3\\times 10^{-10}$)",
        "58.5% sign ($r=+0.42, P=1.7\\times 10^{-6}$)",
        "64.9% sign ($24/37, r=+0.47, P=0.0031$)",
        "57.0% sign ($110/193, r=+0.14, P=0.047$)",
        "66.9% sign ($83/124, P_{\\mathrm{bin}}=2.0\\times 10^{-4}$)",
        "71.0% sign ($44/62, P_{\\mathrm{bin}}=0.0013$)",
    ]
    c_colors = ["#7f8c8d", "#95a5a6", "#1f4e78", "#2980b9", "#d35400", "#2e8b57", "#8e44ad", "#6c3483"]

    yc = np.arange(len(c_labels))
    bars = axC.barh(yc, c_rvals, height=0.62, color=c_colors, alpha=0.92, edgecolor="black", linewidth=0.6)
    axC.axvline(0.0, color="black", linewidth=1.1)
    axC.axhline(1.5, color="#cccccc", linestyle="--", linewidth=0.9)
    axC.axhline(5.5, color="#cccccc", linestyle="--", linewidth=0.9)

    for i, (val, txt) in enumerate(zip(c_rvals, c_sign_txt)):
        if val >= 0:
            axC.text(val + 0.018, i, txt, va="center", ha="left", fontsize=7.5, fontweight="bold", color="#1a1a1a")
        else:
            axC.text(0.018, i, txt, va="center", ha="left", fontsize=7.5, fontweight="bold", color="#555555")

    axC.set_yticks(yc)
    axC.set_yticklabels(c_labels, fontsize=7.9, fontweight="bold")
    axC.invert_yaxis()
    axC.set_xlim(-0.34, 1.05)
    axC.set_xlabel("Signed Effect Correlation ($r$ of $Z_{\\mathrm{TWAS}}$ or Block Allelic $Z$)", fontweight="bold")
    axC.set_title("C. Signed Quantized Block & TWAS Concordance (Options 2A & 2B)", fontweight="bold", loc="left", pad=8)

    # =========================================================================
    # PANEL D: Gene Program Heritability Partition & Viral Variant In Vivo Origins
    # =========================================================================
    axD = fig.add_subplot(gs[1, 1])
    prog_rows = []
    with open(os.path.join(DATA_DIR, "biobank_lcl_gene_program_partition.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            prog_rows.append(r)

    prog_order = [
        ("1_MHC_and_Antigen_Processing_HLA_ERAP", "1. MHC & ERAP\nAntigen Pres.\n[256 eQTLs]"),
        ("2_Biobank_T_NK_Cell_Extrinsic_Immune_Loci", "2. Non-HLA\nT/NK Immune\n[24 eQTLs]"),
        ("4_BCell_Identity_and_Zp_Rp_Lytic_Switch_TFs", "3. B-Cell & Zp\nLytic TFs\n[17 eQTLs]"),
        ("3_Epigenetic_CpG_Methylation_KZNF_Episome_Silencers", "4. Epigenetic\nCpG / KZNF\n[380 eQTLs]"),
        ("6_MAGE_BGLF5_Lytic_Shutoff_Coupled_Program_Top500", "5. MAGE BGLF5\nLytic Top 500\n[341 eQTLs]"),
        ("5_MAGE_LMP1_NFkB_Activated_Host_Program_Top500", "6. MAGE LMP1\nNF-$\\kappa$B Top 500\n[381 eQTLs]"),
    ]

    def get_prog_enrich(prog_key, trait_key):
        for r in prog_rows:
            if r["gene_program"] == prog_key and r["trait"] == trait_key:
                return min(6.0, float(r["fold_enrichment_vs_bg"])), float(r["fold_enrichment_vs_bg"])
        return 1.0, 1.0

    xp = np.arange(len(prog_order))
    bw = 0.20
    tr_specs = [
        ("Biobank_Blood_UKB_AoU", "Biobank Blood ($N\\approx 900\\mathrm{k}$)", "#1f4e78", -1.5 * bw),
        ("Biobank_Saliva_AoU", "Biobank Saliva ($N\\approx 245\\mathrm{k}$)", "#2e8b57", -0.5 * bw),
        ("LCL_RNA_ProdFork_ViralCond", "LCL `RNA_ProdFork_ViralCond`", "#e67e22", +0.5 * bw),
        ("LCL_RNA_LMP1_ViralCond", "LCL `RNA_LMP1_ViralCond`", "#c0392b", +1.5 * bw),
    ]

    for tr_key, tr_lbl, col, offset in tr_specs:
        capped_vals, raw_vals = zip(*[get_prog_enrich(pk, tr_key) for pk, _ in prog_order])
        axD.bar(xp + offset, capped_vals, width=bw, color=col, alpha=0.90, label=tr_lbl, edgecolor="black", linewidth=0.5)
        for idx, (cv, rv) in enumerate(zip(capped_vals, raw_vals)):
            if rv > 6.0:
                axD.text(xp[idx] + offset, 6.08, f"{rv:.0f}x", ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=col)

    axD.axhline(1.0, color="#333333", linestyle="--", linewidth=1.0, label="Non-MHC Background ($1.0\\times$)")

    # Add inset box summarizing Analysis 4B: In Vivo Biobank (UKB+AoU) vs 1kGP Viral Variant Origins
    inset_txt = (
        "In Vivo Biobank (UKB+AoU) vs. 1kGP Viral Variants (112 Loci):\n"
        "• Wild-Type Polymorphic In Vivo & LCLs (12 loci, $\\mathrm{AF}\\geq 5\\%$):\n"
        "  - `12,129:T>A` (IR1/Wp #1 CN QTL): 80.0% UKB, 86.0% AoU\n"
        "  - `36,480:T>C` & `36,494:G>C` (EBNA-2): 32.6% UKB, 44.9% AoU\n"
        "• B95-8 Lab Culture-Specific Mutations (<1.8% in UKB/AoU):\n"
        "  - `65,739:AT>A` (BORF2), `154,091:A>G` (BALF5), `163,858:G>T`,\n"
        "    `91,192:G>A` (BRLF1), `132,965:G>GC`: 0.0%–1.75% in vivo"
    )
    axD.text(
        0.99, 0.78, inset_txt,
        transform=axD.transAxes,
        ha="right", va="top",
        fontsize=6.6,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#f8f9f9", edgecolor="#2c3e50", alpha=0.96),
    )

    axD.set_xticks(xp)
    axD.set_xticklabels([lbl for _, lbl in prog_order], fontsize=7.0, fontweight="bold")
    axD.set_ylabel("Mean $\\chi^2$ Fold-Enrichment vs.\nNon-MHC Expressed Background (capped at 6x)", fontweight="bold")
    axD.set_ylim(0, 7.35)
    axD.set_title("D. Partitioned Program Enrichment & In Vivo Viral Variant Origins", fontweight="bold", loc="left", pad=8)
    axD.legend(loc="upper right", ncol=2, frameon=True, facecolor="white", edgecolor="#cccccc", framealpha=0.95, fontsize=6.9)

    fig.suptitle(
        "In Vivo Biobank EBV GWAS (UKB + All of Us) vs. B-Cell-Intrinsic Viral-Conditioned LCL Architecture",
        fontsize=13.5, fontweight="bold", y=0.98
    )

    png_path = os.path.join(FIG_DIR, "fig8_biobank_vs_lcl_architecture.png")
    pdf_path = os.path.join(FIG_DIR, "fig8_biobank_vs_lcl_architecture.pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()
