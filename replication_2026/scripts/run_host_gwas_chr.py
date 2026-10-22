#!/usr/bin/env python3
"""
Per-chromosome host genetic association runner (Layer B1, Layer B2, and Layer C).
Runs PLINK2 --glm across all 11 EBV WGS/RNA phenotypes (unconditioned and viral-conditioned)
for a single autosome (chr1..chr22) and extracts:
  1. Layer C: Genome-wide common variant (MAF >= 5%) summary stats, exact chi-square
     distribution for lambda_GC / QQ plots, 500-kb Manhattan bins, and all P < 1e-4 hits.
  2. Layer B1: Exact association stats for all MAGE fine-mapped cis-eQTL lead variants
     on this chromosome (16,008 genome-wide lead eQTLs).
  3. Layer B2: Window-based association summary (+/- 100 kb) for 55 canonical host
     EBV-regulator genes (entry receptors, HLA class I/II, Latency III B-cell TFs,
     Zp/Rp lytic switch TFs, SMC5/6 & HUSH/PML episome restriction, and DNA methyltransferases).
"""
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
import numpy as np

BASE = "/work/nfp8/EBV"
GWAS_DIR = f"{BASE}/host_gwas"
PLINK2 = "/hpc/group/parrishlab/ms1299/miniconda3/envs/py311_bash5_2/bin/plink2"

# 55 Canonical Host Regulators of EBV Entry, Latency III, Lytic Switch, & Episome Restriction (GRCh38)
CANDIDATE_55_GENES = [
    # Viral entry & B-cell surface co-receptors
    ("CR2", "1", 207454206, 207489895, "Entry_Receptor"),
    ("CD19", "16", 28931971, 28939342, "Entry_CoReceptor"),
    ("CD81", "11", 2376177, 2397419, "Entry_CoReceptor"),
    ("CD40", "20", 46118243, 46129863, "LMP1_Mimic_Receptor"),
    # MHC / HLA class I & II antigen presentation & EBV entry co-factor (HLA class II binds gp42)
    ("HLA-A", "6", 29941260, 29945884, "MHC_HLA_ClassI"),
    ("HLA-B", "6", 31353872, 31357187, "MHC_HLA_ClassI"),
    ("HLA-C", "6", 31268749, 31272130, "MHC_HLA_ClassI"),
    ("HLA-DRA", "6", 32439878, 32445046, "MHC_HLA_ClassII_gp42"),
    ("HLA-DRB1", "6", 32578769, 32589848, "MHC_HLA_ClassII_gp42"),
    ("HLA-DQA1", "6", 32637406, 32654846, "MHC_HLA_ClassII_gp42"),
    ("HLA-DQB1", "6", 32659467, 32668383, "MHC_HLA_ClassII_gp42"),
    ("HLA-DPB1", "6", 33075990, 33089696, "MHC_HLA_ClassII_gp42"),
    ("CIITA", "16", 10877363, 10925178, "MHC_ClassII_Master_TF"),
    # Latency III (Cp / Wp / LMP1p) B-cell identity & NF-kB TFs
    ("RBPJ", "4", 26163459, 26434880, "EBNA2_Cp_LMP1p_Anchor_TF"),
    ("PAX5", "9", 36833272, 37034476, "Wp_Latency_Bcell_TF"),
    ("EBF1", "5", 158695916, 159099923, "EBNA2_Pioneer_Bcell_TF"),
    ("SPI1", "11", 47354860, 47395640, "LMP1p_PU1_TF"),
    ("POU2F2", "19", 42086399, 42139874, "OCT2_Lytic_Repressor_TF"),
    ("POU2AF1", "11", 111351478, 111383439, "OCAB_BOB1_Bcell_Coactivator"),
    ("IRF4", "6", 391752, 411447, "EBNA3C_LMP1_Bcell_Survival_TF"),
    ("IRF8", "16", 85899165, 85922606, "Bcell_Lineage_TF"),
    ("NFKB1", "4", 102501329, 102617302, "LMP1_Canonical_NFkB"),
    ("NFKB2", "10", 102394110, 102402529, "LMP1_NonCanonical_NFkB"),
    ("REL", "2", 60881492, 60931612, "LMP1_NFkB_cRel"),
    ("RELA", "11", 65653596, 65662972, "LMP1_NFkB_p65"),
    ("RELB", "19", 45001449, 45038198, "LMP1_NonCanonical_NFkB"),
    ("TRAF1", "9", 120902397, 120929244, "LMP1_CTAR1_Signaling"),
    ("STAT3", "17", 42313324, 42388442, "Lytic_Refractoriness_TF"),
    ("MYC", "8", 127735434, 127742951, "Lytic_Switch_Repressor"),
    # Lytic Switch (Zp / Rp) Plasma-Cell Differentiation & Stress TFs
    ("PRDM1", "6", 106086320, 106109939, "BLIMP1_Zp_Lytic_Activator"),
    ("XBP1", "22", 28794555, 28800597, "XBP1s_Zp_Lytic_Activator"),
    ("MEF2D", "1", 156461591, 156498761, "Zp_MEF2_Chromatin_Switch"),
    ("KLF4", "9", 107484852, 107490482, "Zp_Rp_Lytic_Activator"),
    ("SP1", "12", 53380176, 53416446, "Zp_Promoter_TF"),
    ("CREB1", "2", 207529892, 207605989, "Zp_ZII_Element_TF"),
    ("ATF1", "12", 50763901, 50821275, "Zp_ZII_Element_TF"),
    ("JUN", "1", 58780791, 58784047, "Zp_AP1_Element_TF"),
    ("FOS", "14", 75278826, 75282230, "Zp_AP1_Element_TF"),
    ("HIF1A", "14", 61695513, 61748259, "Zp_Hypoxia_Lytic_Switch"),
    # Episome Maintenance, Cohesin/CTCF Loop Architecture, & SMC5/6 / PML / HUSH Restriction
    ("CTCF", "16", 67562407, 67639185, "Episome_Chromatin_Insulator"),
    ("PARP1", "1", 226360210, 226408154, "OriP_EBNA1_Modifier"),
    ("SMC5", "9", 70259516, 70345558, "SMC5_6_Episome_Restriction"),
    ("SMC6", "2", 17700110, 17786714, "SMC5_6_Episome_Restriction"),
    ("NSMCE2", "8", 125090808, 125164207, "SMC5_6_SUMO_Ligase"),
    ("PML", "15", 73994673, 74047827, "ND10_Promyelocytic_Body"),
    ("SP100", "2", 230417873, 230552645, "ND10_Episome_Repressor"),
    ("DAXX", "6", 33318558, 33323016, "BNRF1_Target_H33_Chaperone"),
    ("ATRX", "X", 77504880, 77786216, "BNRF1_Target_Chromatin_Remodeler"),
    ("TASOR", "3", 52579368, 52627391, "HUSH_Complex_Episome_Silencer"),
    ("MORC2", "22", 30923476, 30979263, "HUSH_Complex_Effector"),
    # Host Epigenetic CpG Methylation Machinery (Required for BZLF1 Methylated-DNA Binding)
    ("DNMT1", "19", 10133345, 10231286, "Episome_CpG_Maintenance_Methylase"),
    ("DNMT3A", "2", 25227855, 25342590, "Episome_DeNovo_Methylase"),
    ("DNMT3B", "20", 32762385, 32809356, "Episome_DeNovo_Methylase"),
    ("TET2", "4", 105145875, 105279803, "Episome_5mC_Demethylase"),
    ("TP53", "17", 7661779, 7687550, "Host_Cell_Cycle_Checkpoint"),
]

PHENO_NAMES = [
    "WGS_CN_PopAdj", "WGS_CN_ViralCond", "WGS_IR1_PopAdj",
    "RNA_BGLF5_ViralCond", "RNA_BZLF1_ViralCond", "RNA_IELatRatio_ViralCond",
    "RNA_ProdFork_ViralCond", "RNA_LMP1_ViralCond", "RNA_Latent_ViralCond",
    "RNA_BGLF5_Uncond", "RNA_LMP1_Uncond"
]


def main():
    chrom = str(sys.argv[1])
    out_dir = f"{GWAS_DIR}/chr_results"
    os.makedirs(out_dir, exist_ok=True)

    # Check if already done
    done_flag = f"{out_dir}/chr{chrom}.done"
    if os.path.exists(done_flag):
        print(f"chr{chrom} already completed.")
        return

    tmp_dir = f"/tmp/nfp8_gwas_chr{chrom}_{os.getpid()}"
    os.makedirs(tmp_dir, exist_ok=True)
    vcf_path = f"{tmp_dir}/chr{chrom}.vcf.gz"
    url = f"https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

    print(f"[chr{chrom}] Downloading VCF from UCSC: {url}")
    subprocess.run(["curl", "-sSL", "--retry", "10", "--retry-connrefused", "--retry-delay", "5", "-o", vcf_path, url], check=True)

    print(f"[chr{chrom}] Running PLINK2 --glm across 11 phenotypes...")
    cmd = [
        PLINK2,
        "--vcf", vcf_path,
        "--vcf-half-call", "missing",
        "--pheno", f"{GWAS_DIR}/ebv_host_phenotypes.plink.tsv",
        "--no-input-missing-phenotype",
        "--maf", "0.01",
        "--max-alleles", "2",
        "--set-all-var-ids", "@:#:$r:$a",
        "--new-id-max-allele-len", "100", "truncate",
        "--rm-dup", "force-first",
        "--glm", "allow-no-covars", "omit-ref",
        "--threads", "4",
        "--out", f"{tmp_dir}/res"
    ]
    subprocess.run(cmd, check=True)

    # Remove downloaded VCF immediately to free /tmp space
    if os.path.exists(vcf_path):
        os.remove(vcf_path)

    # Load MAGE fine-mapped cis-eQTL lead variants on this chromosome
    mage_eqtls = defaultdict(list)
    with open(f"{GWAS_DIR}/mage_finemapped_lead_eqtls.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["variantChrom"] == f"chr{chrom}" or row["variantChrom"] == chrom:
                mage_eqtls[row["variant_kgpID"]].append(row)

    # Filter candidate 55 genes on this chromosome (+/- 100 kb window)
    chr_cands = [g for g in CANDIDATE_55_GENES if g[1] == chrom]

    for pheno in PHENO_NAMES:
        glm_file = f"{tmp_dir}/res.{pheno}.glm.linear"
        if not os.path.exists(glm_file):
            continue

        n_maf01 = 0
        n_common = 0
        n_gw_5e8 = 0
        n_sug_1e6 = 0
        n_sug_1e5 = 0
        chisq_hist = np.zeros(10000, dtype=np.int64)  # bins of 0.005 on [0, 50]
        high_chisq = []  # exact chi2 >= 15.0 for upper QQ tail
        min_p_row = None
        min_p_val = 2.0

        suggestive_rows = []
        manhattan_bins = {}  # bin_idx (500kb) -> (p_val, row_tuple)
        eqtl_hits = []
        cand_best = {g[0]: {"n_maf01": 0, "n_maf05": 0, "min_p": 2.0, "best_row": None} for g in chr_cands}

        with open(glm_file) as f:
            hdr = f.readline().strip().split("\t")
            # Standard PLINK2 glm.linear columns:
            # #CHROM POS ID REF ALT PROVISIONAL_REF? A1 OMITTED A1_FREQ TEST OBS_CT BETA SE T_STAT P ERRCODE
            col_idx = {h.lstrip("#"): i for i, h in enumerate(hdr)}
            i_pos = col_idx["POS"]
            i_id = col_idx["ID"]
            i_a1 = col_idx["A1"]
            i_freq = col_idx["A1_FREQ"]
            i_n = col_idx["OBS_CT"]
            i_beta = col_idx["BETA"]
            i_se = col_idx["SE"]
            i_t = col_idx["T_STAT"]
            i_p = col_idx["P"]

            for line in f:
                p = line.rstrip("\n").split("\t")
                p_str = p[i_p]
                if p_str == "NA":
                    continue
                try:
                    pval = float(p_str)
                    tstat = float(p[i_t])
                    freq = float(p[i_freq])
                    pos = int(p[i_pos])
                except ValueError:
                    continue

                vid = p[i_id]
                is_common = 0.05 <= freq <= 0.95
                n_maf01 += 1

                if is_common:
                    n_common += 1
                    chi2 = tstat * tstat
                    b_idx = int(chi2 * 200.0)
                    if 0 <= b_idx < 10000:
                        chisq_hist[b_idx] += 1
                    if chi2 >= 15.0:
                        high_chisq.append(round(chi2, 4))
                    if pval < 5e-8:
                        n_gw_5e8 += 1
                    if pval < 1e-6:
                        n_sug_1e6 += 1
                    if pval < 1e-5:
                        n_sug_1e5 += 1
                    if pval < min_p_val:
                        min_p_val = pval
                        min_p_row = [chrom, str(pos), vid, p[i_a1], f"{freq:.4f}", p[i_n], p[i_beta], p[i_se], p[i_t], p_str]

                    m_bin = pos // 500000
                    if m_bin not in manhattan_bins or pval < manhattan_bins[m_bin][0]:
                        manhattan_bins[m_bin] = (pval, [chrom, str(pos), vid, f"{freq:.4f}", p[i_beta], p[i_t], p_str])

                if pval < 1e-4:
                    suggestive_rows.append([chrom, str(pos), vid, p[i_a1], f"{freq:.4f}", p[i_n], p[i_beta], p[i_se], p[i_t], p_str])

                if vid in mage_eqtls:
                    for eq in mage_eqtls[vid]:
                        eqtl_hits.append([
                            chrom, str(pos), vid, eq["variant_rsID"], eq["geneSymbol"],
                            eq["ensemblID"], f"{float(eq['variantPIP']):.4f}", f"{float(eq['log2_aFC']):.4f}",
                            p[i_a1], f"{freq:.4f}", p[i_n], p[i_beta], p[i_se], p[i_t], p_str
                        ])

                for g_name, _, g_lo, g_hi, g_cat in chr_cands:
                    if (g_lo - 100000) <= pos <= (g_hi + 100000):
                        cb = cand_best[g_name]
                        cb["n_maf01"] += 1
                        if is_common:
                            cb["n_maf05"] += 1
                        if pval < cb["min_p"]:
                            cb["min_p"] = pval
                            cb["best_row"] = [
                                g_name, g_cat, chrom, str(pos), vid, p[i_a1],
                                f"{freq:.4f}", p[i_n], p[i_beta], p[i_se], p[i_t], p_str
                            ]

        # Save outputs for this chromosome & phenotype
        summary_obj = {
            "chrom": chrom,
            "pheno": pheno,
            "n_maf01": n_maf01,
            "n_common_maf05": n_common,
            "n_gw_5e8": n_gw_5e8,
            "n_sug_1e6": n_sug_1e6,
            "n_sug_1e5": n_sug_1e5,
            "min_p_row": min_p_row,
            "chisq_hist": chisq_hist.tolist(),
            "high_chisq": high_chisq,
        }
        with open(f"{out_dir}/chr{chrom}_{pheno}_summary.json", "w") as out:
            json.dump(summary_obj, out)

        with open(f"{out_dir}/chr{chrom}_{pheno}_suggestive.tsv", "w") as out:
            out.write("chrom\tpos\tid\ta1\ta1_freq\tn\tbeta\tse\tt_stat\tp\n")
            for r in suggestive_rows:
                out.write("\t".join(r) + "\n")

        with open(f"{out_dir}/chr{chrom}_{pheno}_manhattan.tsv", "w") as out:
            out.write("chrom\tpos\tid\ta1_freq\tbeta\tt_stat\tp\n")
            for b_idx in sorted(manhattan_bins):
                out.write("\t".join(manhattan_bins[b_idx][1]) + "\n")

        with open(f"{out_dir}/chr{chrom}_{pheno}_mage_eqtls.tsv", "w") as out:
            out.write("chrom\tpos\tid\trsid\tgene\tensembl_id\tpip\tlog2_afc\ta1\ta1_freq\tn\tbeta\tse\tt_stat\tp\n")
            for r in eqtl_hits:
                out.write("\t".join(r) + "\n")

        with open(f"{out_dir}/chr{chrom}_{pheno}_candidate55.tsv", "w") as out:
            out.write("gene\tcategory\tn_maf01\tn_maf05\tchrom\tpos\tid\ta1\ta1_freq\tn\tbeta\tse\tt_stat\tp\n")
            for g_name, _, _, _, g_cat in chr_cands:
                cb = cand_best[g_name]
                if cb["best_row"] is not None:
                    br = cb["best_row"]
                    out.write("\t".join([br[0], br[1], str(cb["n_maf01"]), str(cb["n_maf05"])] + br[2:]) + "\n")

        os.remove(glm_file)

    shutil.rmtree(tmp_dir, ignore_errors=True)
    with open(done_flag, "w") as out:
        out.write("OK\n")
    print(f"[chr{chrom}] Completed all 11 phenotypes successfully.")


if __name__ == "__main__":
    main()
