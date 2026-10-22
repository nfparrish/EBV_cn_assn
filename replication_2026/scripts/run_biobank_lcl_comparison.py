#!/usr/bin/env python3
"""
Comprehensive Statistical Comparison of In Vivo Biobank EBV GWAS
(Nyeo/Lareau et al. 2026, Loh/Kamitaki et al. Nature 2026, Schmidt/Okada et al. Nature 2026)
Against Viral-Conditioned B-Cell-Intrinsic LCL Host Genetics (1kGP N=2,499 WGS + MAGE N=702 RNA).

Executes all 4 analyses:
  1. Unsigned Genome-Wide Quantized Block Test (Full Genome vs. Immune/MHC-Masked Sub-Threshold Tail)
     across ~5,457 contiguous 500-kb autosomal blocks and 3,650 cis-eQTL blocks.
  2. Signed Quantized Block Tests (Both Options):
     - Option 2A: Allelic-Sign Concordance & Precision-Weighted Z-Regression across Quantized
       Biobank Blocks (both Biobank-selected block lead variants and LCL-selected block lead variants,
       before and after masking MHC + genome-wide significant immune peaks).
     - Option 2B: Germline TWAS / cis-eQTL Signed Concordance across 12,544 fine-mapped MAGE
       cis-eQTLs (Z_TWAS = sign(log2_aFC) * Z_GWAS), testing cross-study TWAS agreement and
       agreement with empirical MAGE host-viral mRNA coupling (sign(r_RNA)).
  3. Partitioning Sub-Threshold Heritability (Mean Chi-Square & Enrichment) Across 7 Host & Viral
     Gene Expression Programs (MHC/Antigen Processing, Non-HLA Biobank Immune/T-NK Loci,
     Epigenetic CpG Methylation & Heterochromatin Silencers, B-Cell/Lytic Switch TFs,
     MAGE LMP1/NF-kB Targets, MAGE BGLF5 Host-Shutoff Targets, and Background Genes).
  4. Biobank High-Prior Lookup Table (exact Biobank lead SNP via remote tabix + +/-100kb window)
     & In Vivo Biobank (UKB + AoU) vs. 1kGP LCL Viral Variant Architecture Comparison.
"""
import csv
import gzip
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
from collections import defaultdict
import numpy as np
from scipy import stats

BASE = "/work/nfp8/EBV"
GWAS_DIR = f"{BASE}/host_gwas"
CHR_DIR = f"{GWAS_DIR}/chr_results"
BIOBANK_DIR = f"{GWAS_DIR}/biobank"
LAREAU_DIR = f"{GWAS_DIR}/clareaulab_ebv_biobank_gwas"
OUT_DIR = f"{GWAS_DIR}/summary"
TABIX = "/hpc/group/parrishlab/ms1299/miniconda3/envs/py311_bash5_2/bin/tabix"

BLOOD_FILE = f"{BIOBANK_DIR}/EBV_invnorm_UKB_AoU_SE1.CHR_BP.tbl.gz"
SALIVA_FILE = f"{BIOBANK_DIR}/EBV_masked_invnorm.bgen.TOPMed_MAF_0.001.AC_40.stats.gz"

LCL_PHENOS = [
    "WGS_CN_ViralCond",
    "WGS_CN_PopAdj",
    "RNA_IELatRatio_ViralCond",
    "RNA_ProdFork_ViralCond",
    "RNA_BGLF5_ViralCond",
    "RNA_BZLF1_ViralCond",
    "RNA_LMP1_ViralCond",
    "RNA_Latent_ViralCond",
]

# High-prior Biobank & Canonical EBV Host Loci (GRCh38 coordinates)
PRIORITY_LOCI = [
    # MHC Class I & II Antigen Presentation & Entry (chr6:25-35 Mb)
    ("HLA-A", "6", 29941260, 29945884, "MHC_ClassI_Antigen_Presentation"),
    ("HLA-C", "6", 31268749, 31272130, "MHC_ClassI_Antigen_Presentation"),
    ("HLA-B", "6", 31353872, 31357187, "MHC_ClassI_Antigen_Presentation"),
    ("HLA-DRA", "6", 32439878, 32445046, "MHC_ClassII_gp42_CoReceptor"),
    ("HLA-DRB1", "6", 32578769, 32589848, "MHC_ClassII_gp42_CoReceptor"),
    ("HLA-DQA1", "6", 32637406, 32654846, "MHC_ClassII_gp42_CoReceptor"),
    ("HLA-DQB1", "6", 32659467, 32668383, "MHC_ClassII_gp42_CoReceptor"),
    ("HLA-DPB1", "6", 33075990, 33089696, "MHC_ClassII_gp42_CoReceptor"),
    # Non-HLA Antigen Processing Peptidases (Top Biobank Non-MHC Hits)
    ("ERAP1", "5", 96760810, 96808106, "ER_Aminopeptidase_MHC_I_Trimming"),
    ("ERAP2", "5", 96875939, 96919716, "ER_Aminopeptidase_MHC_I_Trimming"),
    ("LNPEP", "5", 96934866, 97037530, "Endosomal_Peptidase_CrossPresentation"),
    ("CIITA", "16", 10877363, 10925178, "MHC_ClassII_Master_Transactivator"),
    # Cell-Extrinsic T/NK Immune Surveillance & Co-Stimulation (Biobank GWAS/ExWAS Hits)
    ("PTPN22", "1", 113813811, 113871761, "TCR_Signaling_Phosphatase_T_NK"),
    ("LY9", "1", 160795836, 160830445, "SLAM_Family_T_NK_B_CrossTalk"),
    ("CTLA4", "2", 203867771, 203873965, "T_Cell_Inhibitory_Checkpoint"),
    ("ICOS", "2", 203936748, 203961577, "T_Cell_CoStimulatory_Receptor"),
    ("NBEAL2", "3", 46979642, 47009695, "NK_CTL_Granule_Exocytosis"),
    ("SLC22A4", "5", 132294312, 132344068, "Ergothioneine_Transporter_Immune_Locus"),
    ("BACH2", "6", 89926528, 90296995, "T_B_Cell_Memory_Differentiation_TF"),
    ("CD27", "12", 6442152, 6451728, "T_NK_CoStimulatory_Receptor_for_CD70"),
    ("IFNG", "12", 68154768, 68159740, "Antiviral_Type_II_Interferon_T_NK"),
    ("SH2B3", "12", 111405926, 111451624, "Cytokine_JAK_STAT_Negative_Regulator"),
    ("ITGAL", "16", 30472658, 30523185, "LFA1_Integrin_Immunological_Synapse"),
    ("GSDMB", "17", 39904595, 39919854, "Gasdermin_B_Granzyme_A_Pyroptosis"),
    ("CD226", "18", 69863420, 69962158, "DNAM1_NK_CD8_Activating_Receptor"),
    ("CD70", "19", 6583183, 6591003, "B_Cell_Ligand_for_CD27_CTL_Surveillance"),
    ("FUT2", "19", 48695965, 48705946, "Mucosal_Secretor_Status_Saliva_Shedding"),
    # B-Cell Intrinsic Entry, Latency III & Lytic Switch TFs
    ("CR2", "1", 207454206, 207489895, "EBV_gp350_Entry_Receptor_CD21"),
    ("EBF1", "5", 158695916, 159099923, "B_Cell_Pioneer_TF_EBNA2_CoFactor"),
    ("IRF4", "6", 391752, 411447, "B_Cell_Survival_TF_EBNA3C_LMP1"),
    ("PRDM1", "6", 106086320, 106109939, "BLIMP1_Plasma_Cell_Zp_Activator"),
    ("PAX5", "9", 36833272, 37034476, "B_Cell_Identity_TF_Wp_Binding_Factor"),
    ("STAT3", "17", 42313324, 42388442, "Lytic_Refractoriness_TF"),
    ("CD40", "20", 46118243, 46129863, "B_Cell_CoReceptor_Mimicked_by_LMP1"),
    # B-Cell Intrinsic Epigenetic CpG Methylation & Heterochromatin Silencers (Unmasked in LCLs)
    ("SMC6", "2", 17700110, 17786714, "SMC5_6_Episome_Restriction_Complex"),
    ("TET2", "4", 105145875, 105279803, "Episome_5mC_Dioxygenase"),
    ("SMC5", "9", 70259516, 70345558, "SMC5_6_Episome_Restriction_Complex"),
    ("MED21", "12", 27016900, 27029600, "Mediator_Subunit_LMP1_Transactivation"),
    ("DNMT1", "19", 10133345, 10231286, "Maintenance_CpG_Methyltransferase"),
    ("ZNF814", "19", 57862000, 57891000, "KRAB_ZNF_Heterochromatin_Silencer"),
    ("DNMT3B", "20", 32762385, 32809356, "DeNovo_CpG_Methyltransferase_ZRE_Switch"),
]


def load_lcl_targets():
    """Load all 12,544 MAGE cis-eQTLs and all 5,457 LCL 500-kb Manhattan block leads."""
    # Map (chrom, pos) -> set of (ref, alt) needed for exact matching
    needed_positions = defaultdict(set)

    # 1. Load MAGE cis-eQTLs across all 8 primary LCL phenotypes
    eqtl_meta = {}  # vid -> dict of meta
    eqtl_lcl = defaultdict(dict)  # vid -> {pheno: (t_stat, p_val, beta, se)}
    for c in range(1, 23):
        chrom = str(c)
        for ph in LCL_PHENOS:
            fpath = f"{CHR_DIR}/chr{chrom}_{ph}_mage_eqtls.tsv"
            if not os.path.exists(fpath):
                continue
            with open(fpath) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    vid = row["id"]
                    parts = vid.split(":")
                    pos = int(parts[1])
                    ref, alt = parts[2], parts[3]
                    if vid not in eqtl_meta:
                        eqtl_meta[vid] = {
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "rsid": row["rsid"],
                            "gene": row["gene"],
                            "ensembl_id": row["ensembl_id"].split(".")[0],
                            "pip": float(row["pip"]),
                            "log2_afc": float(row["log2_afc"]),
                            "a1_freq": float(row["a1_freq"]),
                        }
                        needed_positions[(chrom, pos)].add((ref, alt, vid))
                    # Verify only highest PIP per (vid, gene) or store by vid
                    eqtl_lcl[vid][ph] = (
                        float(row["t_stat"]),
                        float(row["p"]),
                        float(row["beta"]),
                        float(row["se"]),
                    )

    #Also store gene-specific eQTL list because one variant can be lead eQTL for >1 gene
    eqtl_gene_rows = []
    seen_vg = set()
    for c in range(1, 23):
        chrom = str(c)
        fpath = f"{CHR_DIR}/chr{chrom}_WGS_CN_ViralCond_mage_eqtls.tsv"
        with open(fpath) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                vid = row["id"]
                gene = row["gene"]
                if (vid, gene) in seen_vg:
                    continue
                seen_vg.add((vid, gene))
                parts = vid.split(":")
                eqtl_gene_rows.append({
                    "vid": vid,
                    "chrom": chrom,
                    "pos": int(parts[1]),
                    "ref": parts[2],
                    "alt": parts[3],
                    "rsid": row["rsid"],
                    "gene": gene,
                    "ensembl_id": row["ensembl_id"].split(".")[0],
                    "pip": float(row["pip"]),
                    "log2_afc": float(row["log2_afc"]),
                    "a1_freq": float(row["a1_freq"]),
                })

    # 2. Load 500-kb Manhattan block leads across all 8 primary LCL phenotypes
    # Key: (chrom, bin_idx), where bin_idx = pos // 500000
    lcl_blocks = defaultdict(dict)  # (chrom, bin_idx) -> {pheno: (vid, pos, ref, alt, t_stat, p_val, beta, freq)}
    for c in range(1, 23):
        chrom = str(c)
        for ph in LCL_PHENOS:
            fpath = f"{CHR_DIR}/chr{chrom}_{ph}_manhattan.tsv"
            if not os.path.exists(fpath):
                continue
            with open(fpath) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    vid = row["id"]
                    parts = vid.split(":")
                    pos = int(row["pos"])
                    ref, alt = parts[2], parts[3]
                    b_idx = pos // 500000
                    lcl_blocks[(chrom, b_idx)][ph] = (
                        vid,
                        pos,
                        ref,
                        alt,
                        float(row["t_stat"]),
                        float(row["p"]),
                        float(row["beta"]),
                        float(row["a1_freq"]),
                    )
                    needed_positions[(chrom, pos)].add((ref, alt, vid))

    return needed_positions, eqtl_meta, eqtl_gene_rows, eqtl_lcl, lcl_blocks


def parse_blood_worker(needed_positions):
    """Stream UKB+AoU Blood summary stats and extract:
    1. Exact matches to needed_positions (MAGE eQTLs + LCL block leads).
    2. 500-kb block summaries (max_chisq, sum_chisq, n_snps, n_p05, lead variant).
    3. Lead variant within +/-100 kb of each PRIORITY_LOCI gene.
    4. All genome-wide significant 1-Mb loci (P < 5e-8).
    """
    matched_variants = {}  # vid -> (beta_alt, se, z_alt, pval, direction)
    blocks = {}  # (chrom, b_idx) -> [n_snps, sum_chisq, n_p05, max_chisq, min_p, lead_pos, lead_a1, lead_a2, lead_beta, lead_se, lead_z]
    locus_best = {}  # gene -> (min_p, max_chisq, chrom, pos, a1, a2, beta, se, z)
    gw_peaks = {}  # (chrom, mb_idx) -> (min_p, max_chisq, chrom, pos, a1, a2, beta, se, z)

    # Index priority loci by chromosome for fast lookup
    loci_by_chr = defaultdict(list)
    for g, ch, st, en, cat in PRIORITY_LOCI:
        loci_by_chr[ch].append((g, max(1, st - 100000), en + 100000))

    proc = subprocess.Popen(["gzip", "-dc", BLOOD_FILE], stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
    header = proc.stdout.readline()
    # CHR BP MarkerName Allele1 Allele2 Effect StdErr P-value Direction
    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        chrom = p[0]
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        if not chrom.isdigit():
            continue
        c_int = int(chrom)
        if c_int < 1 or c_int > 22:
            continue
        try:
            pos = int(p[1])
            beta = float(p[5])
            se = float(p[6])
            pval = float(p[7])
        except ValueError:
            continue
        if se <= 0:
            continue
        a1 = p[3].upper()
        a2 = p[4].upper()
        z = beta / se
        chi2 = z * z

        # 1. Check exact needed positions
        targets = needed_positions.get((chrom, pos))
        if targets is not None:
            for ref, alt, vid in targets:
                if a1 == alt and a2 == ref:
                    matched_variants[vid] = (beta, se, z, pval, p[8])
                elif a1 == ref and a2 == alt:
                    matched_variants[vid] = (-beta, se, -z, pval, p[8])

        # 2. Update 500-kb block stats
        b_key = (chrom, pos // 500000)
        b_rec = blocks.get(b_key)
        if b_rec is None:
            blocks[b_key] = [1, chi2, 1 if pval < 0.05 else 0, chi2, pval, pos, a1, a2, beta, se, z]
        else:
            b_rec[0] += 1
            b_rec[1] += chi2
            if pval < 0.05:
                b_rec[2] += 1
            if chi2 > b_rec[3]:
                b_rec[3] = chi2
                b_rec[4] = pval
                b_rec[5] = pos
                b_rec[6] = a1
                b_rec[7] = a2
                b_rec[8] = beta
                b_rec[9] = se
                b_rec[10] = z

        # 3. Check priority loci on this chromosome
        for g, w_st, w_en in loci_by_chr.get(chrom, ()):
            if w_st <= pos <= w_en:
                cur = locus_best.get(g)
                if cur is None or chi2 > cur[1]:
                    locus_best[g] = (pval, chi2, chrom, pos, a1, a2, beta, se, z)

        # 4. Track genome-wide significant peaks (P < 5e-8) in 1-Mb bins
        if pval < 5e-8:
            mb_key = (chrom, pos // 1000000)
            cur = gw_peaks.get(mb_key)
            if cur is None or chi2 > cur[1]:
                gw_peaks[mb_key] = (pval, chi2, chrom, pos, a1, a2, beta, se, z)

    proc.wait()
    return matched_variants, blocks, locus_best, gw_peaks


def parse_saliva_worker(needed_positions):
    """Stream All of Us Saliva EBV summary stats and extract matching structures."""
    matched_variants = {}  # vid -> (beta_alt, se, z_alt, pval, a1freq)
    blocks = {}
    locus_best = {}
    gw_peaks = {}

    loci_by_chr = defaultdict(list)
    for g, ch, st, en, cat in PRIORITY_LOCI:
        loci_by_chr[ch].append((g, max(1, st - 100000), en + 100000))

    proc = subprocess.Popen(["gzip", "-dc", SALIVA_FILE], stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
    header = proc.stdout.readline()
    # SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ INFO CHISQ_LINREG P_LINREG BETA SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM P_BOLT_LMM
    for line in proc.stdout:
        p = line.rstrip("\n").split("\t")
        if len(p) < 16:
            continue
        chrom = p[1]
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        if not chrom.isdigit():
            continue
        c_int = int(chrom)
        if c_int < 1 or c_int > 22:
            continue
        try:
            pos = int(p[2])
            a1freq = float(p[6])
            beta = float(p[10])
            se = float(p[11])
            pval = float(p[15]) if p[15] != "NA" else float(p[13])
            chi2_bolt = float(p[14]) if p[14] != "NA" else float(p[12])
        except ValueError:
            continue
        if se <= 0:
            continue
        a1 = p[4].upper()
        a0 = p[5].upper()
        z_sign = 1.0 if beta >= 0 else -1.0
        z = z_sign * math.sqrt(max(0.0, chi2_bolt))
        chi2 = chi2_bolt

        # 1. Check exact needed positions
        targets = needed_positions.get((chrom, pos))
        if targets is not None:
            for ref, alt, vid in targets:
                if a1 == alt and a0 == ref:
                    matched_variants[vid] = (beta, se, z, pval, a1freq)
                elif a1 == ref and a0 == alt:
                    matched_variants[vid] = (-beta, se, -z, pval, 1.0 - a1freq)

        # 2. Update 500-kb block stats (restrict to common variants MAF >= 0.01 for clean comparison)
        if 0.01 <= a1freq <= 0.99:
            b_key = (chrom, pos // 500000)
            b_rec = blocks.get(b_key)
            if b_rec is None:
                blocks[b_key] = [1, chi2, 1 if pval < 0.05 else 0, chi2, pval, pos, a1, a0, beta, se, z]
            else:
                b_rec[0] += 1
                b_rec[1] += chi2
                if pval < 0.05:
                    b_rec[2] += 1
                if chi2 > b_rec[3]:
                    b_rec[3] = chi2
                    b_rec[4] = pval
                    b_rec[5] = pos
                    b_rec[6] = a1
                    b_rec[7] = a0
                    b_rec[8] = beta
                    b_rec[9] = se
                    b_rec[10] = z

        # 3. Check priority loci
        for g, w_st, w_en in loci_by_chr.get(chrom, ()):
            if w_st <= pos <= w_en:
                cur = locus_best.get(g)
                if cur is None or chi2 > cur[1]:
                    locus_best[g] = (pval, chi2, chrom, pos, a1, a0, beta, se, z)

        # 4. Track genome-wide significant peaks (P < 5e-8) in 1-Mb bins
        if pval < 5e-8:
            mb_key = (chrom, pos // 1000000)
            cur = gw_peaks.get(mb_key)
            if cur is None or chi2 > cur[1]:
                gw_peaks[mb_key] = (pval, chi2, chrom, pos, a1, a0, beta, se, z)

    proc.wait()
    return matched_variants, blocks, locus_best, gw_peaks


def test_exact_snps_in_lcls(target_snps):
    """Given a list of (chrom, pos, a1, a2) from Biobank lead variants, query UCSC 1kGP GRCh38 VCFs
    via tabix and compute exact linear regression against our LCL phenotypes!"""
    # Load LCL phenotypes (already residualized for covariates in ebv_host_phenotypes.plink.tsv!)
    pheno_data = {}  # iid -> {ph: val}
    with open(f"{GWAS_DIR}/ebv_host_phenotypes.plink.tsv") as f:
        hdr = f.readline().strip().split("\t")
        ph_cols = hdr[1:]
        for line in f:
            p = line.strip().split("\t")
            iid = p[0]
            pheno_data[iid] = {}
            for idx, ph in enumerate(ph_cols):
                if idx + 1 < len(p) and p[idx + 1] != "NA":
                    pheno_data[iid][ph] = float(p[idx + 1])

    # Load covariates to perform exact Frisch-Waugh-Lovell residualization on both G and Y
    META_DIR = "/work/nfp8/EBV/metadata"
    ped = {}
    if os.path.exists(f"{META_DIR}/ped3202.txt"):
        with open(f"{META_DIR}/ped3202.txt") as f:
            for line in f:
                if line.startswith("FamilyID"):
                    continue
                p = line.split()
                if len(p) >= 7:
                    ped[p[1]] = {"pop": p[5], "sex": p[4]}
    host_pcs = {}
    if os.path.exists(f"{GWAS_DIR}/1kgp_top_pcs.tsv"):
        with open(f"{GWAS_DIR}/1kgp_top_pcs.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                host_pcs[r["sampleID"]] = [float(r[f"PC{i}"]) for i in range(1, 11)]
    vcov = {}
    if os.path.exists(f"{META_DIR}/all_viral_covars_3200.tsv"):
        with open(f"{META_DIR}/all_viral_covars_3200.tsv") as f:
            hdr_v = f.readline().strip().split("\t")
            for line in f:
                p = line.strip().split("\t")
                vcov[p[0]] = dict(zip(hdr_v, p))
    gcbias = {}
    if os.path.exists(f"{META_DIR}/all_gcbias_2502.tsv"):
        with open(f"{META_DIR}/all_gcbias_2502.tsv") as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 5:
                    try:
                        gcbias[p[0]] = float(p[4])
                    except ValueError:
                        pass
    by_donor_rna = {}
    if os.path.exists(f"{META_DIR}/all_mage_ebv_rna_779.tsv"):
        with open(f"{META_DIR}/all_mage_ebv_rna_779.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                by_donor_rna.setdefault(r["donor"], []).append(r)

    viral_cols = [
        "v12129_TA", "v163858_GT", "v91192_GA", "v132965_GGC",
        "v162476_TC", "v169200_CT", "v38566_GA", "v65013_CT",
        "v81920_CT", "v103742_AG", "v147913_GA"
    ]

    def build_fwl_covars(samples, include_batch=False):
        if not ped or not vcov:
            return np.ones((len(samples), 1), dtype=np.float64)
        n = len(samples)
        cols = [np.ones(n, dtype=np.float64)]
        pops = [ped.get(s, {}).get("pop", "UNK") for s in samples]
        for u in sorted(set(pops))[1:]:
            cols.append(np.array([1.0 if g == u else 0.0 for g in pops], dtype=np.float64))
        sex = [ped.get(s, {}).get("sex", "0") for s in samples]
        for u in sorted(set(sex))[1:]:
            cols.append(np.array([1.0 if g == u else 0.0 for g in sex], dtype=np.float64))
        if include_batch and by_donor_rna:
            batches = [str(by_donor_rna.get(s, [{"batch": "0"}])[0]["batch"]) for s in samples]
            for u in sorted(set(batches))[1:]:
                cols.append(np.array([1.0 if g == u else 0.0 for g in batches], dtype=np.float64))
        for k in range(10):
            cols.append(np.array([host_pcs.get(s, [0.0] * 10)[k] for s in samples], dtype=np.float64))
        ir1 = np.array([
            float(vcov[s]["ir1_copies"]) if vcov.get(s, {}).get("ir1_copies", "NA") != "NA" and float(vcov[s]["ir1_copies"]) > 0 else np.nan
            for s in samples
        ], dtype=np.float64)
        ir1[np.isnan(ir1)] = np.nanmean(ir1) if not np.all(np.isnan(ir1)) else 0.0
        cols.append(ir1)
        gc = np.array([gcbias.get(s, np.nan) for s in samples], dtype=np.float64)
        gc[np.isnan(gc)] = np.nanmean(gc) if not np.all(np.isnan(gc)) else 0.0
        cols.append(gc)
        for vc in viral_cols:
            arr = np.array([
                float(vcov[s][vc]) if vcov.get(s, {}).get(vc, "NA") not in ("NA", "") else np.nan
                for s in samples
            ], dtype=np.float64)
            arr[np.isnan(arr)] = np.nanmean(arr) if not np.all(np.isnan(arr)) else 0.0
            cols.append(arr)
        return np.column_stack(cols)

    import pickle, tempfile
    cache_exact = f"{BIOBANK_DIR}/exact_lcl_fwl_cache.pkl"
    if os.path.exists(cache_exact):
        with open(cache_exact, "rb") as pf:
            return pickle.load(pf)

    by_chr = defaultdict(list)
    for chrom, pos, a1, a2 in target_snps:
        by_chr[str(chrom)].append((int(pos), a1.upper(), a2.upper()))

    results = {}  # (chrom, pos) -> {ph: (beta_a1, se, t_stat, pval, a1_freq, n)}
    if True:
        for chrom, snp_list in sorted(by_chr.items(), key=lambda x: int(x[0])):
            url = f"https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
            regions = [f"{chrom}:{pos}-{pos}" for pos, _, _ in sorted(set(snp_list))]
            cmd = [TABIX, "-h", url] + regions
            try:
                out = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL, cwd="/hpc/home/nfp8", timeout=15)
            except Exception as e:
                print(f"[WARN] tabix failed on chr{chrom}: {e}")
                continue
            samples = []
            wanted_pos = defaultdict(list)
            for pos, a1, a2 in snp_list:
                wanted_pos[pos].append((a1, a2))
            for line in out.splitlines():
                if line.startswith("#CHROM"):
                    samples = line.split("\t")[9:]
                elif not line.startswith("#"):
                    parts = line.split("\t")
                    pos = int(parts[1])
                    ref = parts[3].upper()
                    alt = parts[4].upper()
                    if pos not in wanted_pos:
                        continue
                    matched_alleles = None
                    for a1, a2 in wanted_pos[pos]:
                        if (ref == a1 and alt == a2) or (ref == a2 and alt == a1):
                            matched_alleles = (a1, a2)
                            break
                    if matched_alleles is None:
                        continue
                    a1, a2 = matched_alleles
                    flip = -1.0 if (ref == a1 and alt == a2) else 1.0  # dosage of ALT vs dosage of A1
                    gts = parts[9:]
                    dosages = {}
                    for s, gt in zip(samples, gts):
                        if s not in pheno_data:
                            continue
                        g0 = gt.split(":")[0]
                        if "|" in g0:
                            al = g0.split("|")
                        elif "/" in g0:
                            al = g0.split("/")
                        else:
                            continue
                        if al[0] not in ("0", "1") or al[1] not in ("0", "1"):
                            continue
                        d_alt = int(al[0]) + int(al[1])
                        d_a1 = (2 - d_alt) if flip < 0 else d_alt
                        dosages[s] = float(d_a1)
                    if len(dosages) < 100:
                        continue
                    results[(chrom, pos)] = {}
                    for ph in LCL_PHENOS:
                        s_list, xs, ys = [], [], []
                        for s, d in dosages.items():
                            yv = pheno_data[s].get(ph)
                            if yv is not None:
                                s_list.append(s)
                                xs.append(d)
                                ys.append(yv)
                        if len(xs) < 100:
                            continue
                        x_arr = np.array(xs, dtype=np.float64)
                        y_arr = np.array(ys, dtype=np.float64)
                        af = float(np.mean(x_arr) / 2.0)
                        vx = float(np.var(x_arr, ddof=1))
                        if vx < 1e-6:
                            continue
                        n = len(x_arr)
                        X_cov = build_fwl_covars(s_list, include_batch=ph.startswith("RNA_"))
                        b_x, _, rank_cov, _ = np.linalg.lstsq(X_cov, x_arr, rcond=None)
                        b_y, _, _, _ = np.linalg.lstsq(X_cov, y_arr, rcond=None)
                        xm = x_arr - X_cov @ b_x
                        ym = y_arr - X_cov @ b_y
                        denom_xx = float(np.dot(xm, xm))
                        if denom_xx < 1e-8:
                            continue
                        df = max(1, n - int(rank_cov) - 1)
                        beta = float(np.dot(xm, ym) / denom_xx)
                        resid = ym - beta * xm
                        s2 = float(np.dot(resid, resid) / df)
                        se = math.sqrt(s2 / denom_xx)
                        t_stat = beta / se if se > 0 else 0.0
                        pval = float(2.0 * stats.t.sf(abs(t_stat), df=df))
                        results[(chrom, pos)][ph] = (beta, se, t_stat, pval, af, n)
    with open(cache_exact, "wb") as pf:
        pickle.dump(results, pf)
    return results


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Step 1: Loading LCL eQTLs (12,544) and 500-kb Manhattan block leads (5,457)...")
    needed_positions, eqtl_meta, eqtl_gene_rows, eqtl_lcl, lcl_blocks = load_lcl_targets()
    print(f"  Loaded {len(eqtl_meta)} unique eQTL variants ({len(eqtl_gene_rows)} variant-gene pairs) and {len(lcl_blocks)} 500-kb blocks ({len(needed_positions)} total genomic positions).")

    import pickle
    cache_pkl = f"{BIOBANK_DIR}/parsed_biobank_cache.pkl"
    if os.path.exists(cache_pkl):
        print("Step 2: Loading cached Biobank Blood and Saliva structures from pickle...")
        with open(cache_pkl, "rb") as pf:
            blood_matched, blood_blocks, blood_loci, blood_gw_peaks, saliva_matched, saliva_blocks, saliva_loci, saliva_gw_peaks = pickle.load(pf)
    else:
        print("Step 2: Streaming Biobank Blood (UKB+AoU) and Saliva (AoU) summary stats in parallel...")
        with mp.Pool(2) as pool:
            res_blood_async = pool.apply_async(parse_blood_worker, (needed_positions,))
            res_saliva_async = pool.apply_async(parse_saliva_worker, (needed_positions,))
            blood_matched, blood_blocks, blood_loci, blood_gw_peaks = res_blood_async.get()
            saliva_matched, saliva_blocks, saliva_loci, saliva_gw_peaks = res_saliva_async.get()
        with open(cache_pkl, "wb") as pf:
            pickle.dump((blood_matched, blood_blocks, blood_loci, blood_gw_peaks, saliva_matched, saliva_blocks, saliva_loci, saliva_gw_peaks), pf)

    print(f"  Blood matched {len(blood_matched)} target variants, {len(blood_blocks)} 500-kb blocks, {len(blood_gw_peaks)} 1-Mb GW-sig bins.")
    print(f"  Saliva matched {len(saliva_matched)} target variants, {len(saliva_blocks)} 500-kb blocks, {len(saliva_gw_peaks)} 1-Mb GW-sig bins.")

    # Identify all 500-kb blocks that lie in the MHC (chr6:25-35 Mb) or within +/- 1 Mb of any Biobank P < 5e-8 peak
    immune_gw_peaks = set(blood_gw_peaks.keys()) | set(saliva_gw_peaks.keys())
    blood_gw_set = set(blood_gw_peaks.keys())
    saliva_gw_set = set(saliva_gw_peaks.keys())

    def is_mhc_or_gwsig_block(chrom, b_idx, gw_set=None):
        if gw_set is None:
            gw_set = immune_gw_peaks
        pos_mid = b_idx * 500000 + 250000
        if chrom == "6" and 25000000 <= pos_mid <= 35000000:
            return True
        mb_idx = pos_mid // 1000000
        for dm in (-1, 0, 1):
            if (chrom, mb_idx + dm) in gw_set:
                return True
        return False

    # =========================================================================
    # ANALYSIS 1: UNSIGNED QUANTIZED BLOCK TEST (Unmasked vs. Immune/MHC-Masked)
    # =========================================================================
    print("Step 3: Running Analysis 1 (Unsigned Quantized Block Test: Unmasked vs. Immune/MHC-Masked)...")
    unsigned_rows = []

    for bio_name, bio_blocks, cohort_gw_set in [
        ("Blood_UKB_AoU", blood_blocks, blood_gw_set),
        ("Saliva_AoU", saliva_blocks, saliva_gw_set),
    ]:
        for mask_mode in ["Unmasked_All_Blocks", "Immune_MHC_Masked_SubThreshold"]:
            aligned = []
            for (chrom, b_idx), b_rec in bio_blocks.items():
                if (chrom, b_idx) not in lcl_blocks:
                    continue
                is_masked_any = is_mhc_or_gwsig_block(chrom, b_idx, immune_gw_peaks)
                if mask_mode == "Immune_MHC_Masked_SubThreshold" and is_masked_any:
                    continue
                is_cohort_gwsig_or_mhc = is_mhc_or_gwsig_block(chrom, b_idx, cohort_gw_set)
                n_snps, sum_chi2, n_p05, max_chi2, min_p = b_rec[0], b_rec[1], b_rec[2], b_rec[3], b_rec[4]
                mean_chi2 = sum_chi2 / max(1, n_snps)
                if mask_mode == "Unmasked_All_Blocks" and (is_cohort_gwsig_or_mhc or min_p < 5e-8):
                    tier = "Tier0_GWSig_or_MHC_P_lt_5e-8"
                elif min_p < 1e-5:
                    tier = "Tier1_SubThreshold_5e-8_to_1e-5"
                elif min_p < 1e-3:
                    tier = "Tier2_SubThreshold_1e-5_to_1e-3"
                elif min_p < 0.01:
                    tier = "Tier3_Nominal_1e-3_to_0.01"
                else:
                    tier = "Tier4_Background_P_ge_0.01"
                aligned.append((chrom, b_idx, min_p, max_chi2, mean_chi2, tier, lcl_blocks[(chrom, b_idx)], n_snps))

            for ph in ["WGS_CN_ViralCond", "WGS_CN_PopAdj", "RNA_IELatRatio_ViralCond", "RNA_ProdFork_ViralCond", "RNA_BGLF5_ViralCond", "RNA_LMP1_ViralCond"]:
                pts = [(x[3], x[4], x[2], x[5], x[6][ph][4] ** 2, x[6][ph][5], x[7]) for x in aligned if ph in x[6]]
                bio_max_c = np.array([p[0] for p in pts], dtype=np.float64)
                bio_mean_c = np.array([p[1] for p in pts], dtype=np.float64)
                lcl_max_c = np.array([p[4] for p in pts], dtype=np.float64)
                n_snps_arr = np.array([p[6] for p in pts], dtype=np.float64)
                rho_max, p_rho_max = stats.spearmanr(bio_max_c, lcl_max_c)
                rho_mean, p_rho_mean = stats.spearmanr(bio_mean_c, lcl_max_c)

                # Partial Spearman correlation controlling for block SNP density M_b
                rk_b = stats.rankdata(bio_max_c)
                rk_l = stats.rankdata(lcl_max_c)
                rk_n = stats.rankdata(n_snps_arr)
                rk_b_res = rk_b - np.poly1d(np.polyfit(rk_n, rk_b, 1))(rk_n)
                rk_l_res = rk_l - np.poly1d(np.polyfit(rk_n, rk_l, 1))(rk_n)
                rho_max_adj, p_rho_max_adj = stats.pearsonr(rk_b_res, rk_l_res)

                # Quantile / Decile sign-free concordance using SNP-density-unconfounded bio_mean_c
                med_bio = np.median(bio_mean_c)
                med_lcl = np.median(lcl_max_c)
                concordant_above_below = int(np.sum((bio_mean_c >= med_bio) == (lcl_max_c >= med_lcl)))
                n_tot = len(pts)
                p_binom_med = float(stats.binomtest(concordant_above_below, n_tot, 0.5).pvalue)

                # Top 10% decile overlap using SNP-count-residualized block max chi-square
                q90_bio = np.quantile(rk_b_res, 0.90)
                q90_lcl = np.quantile(rk_l_res, 0.90)
                obs_top10 = int(np.sum((rk_b_res >= q90_bio) & (rk_l_res >= q90_lcl)))
                exp_top10 = n_tot * 0.10 * 0.10
                a = obs_top10
                b = int(np.sum((rk_b_res >= q90_bio) & (rk_l_res < q90_lcl)))
                c = int(np.sum((rk_b_res < q90_bio) & (rk_l_res >= q90_lcl)))
                d = int(np.sum((rk_b_res < q90_bio) & (rk_l_res < q90_lcl)))
                or_top10, p_top10 = stats.fisher_exact([[a, b], [c, d]])

                tier_groups = defaultdict(list)
                for p in pts:
                    tier_groups[p[3]].append(p)
                bg_pool = tier_groups.get("Tier4_Background_P_ge_0.01") or tier_groups.get("Tier3_Nominal_1e-3_to_0.01") or pts
                bg_lcl_mean = float(np.mean([p[4] for p in bg_pool]))

                for t_name in sorted(tier_groups.keys()):
                    t_pts = tier_groups[t_name]
                    t_lcl_c = np.array([p[4] for p in t_pts])
                    t_lcl_p = np.array([p[5] for p in t_pts])
                    t_bio_c = np.array([p[0] for p in t_pts])
                    _, p_mw = stats.mannwhitneyu(t_lcl_c, [p[4] for p in bg_pool], alternative="two-sided") if len(t_pts) >= 5 else (0.0, 1.0)
                    unsigned_rows.append({
                        "biobank_cohort": bio_name,
                        "mask_mode": mask_mode,
                        "lcl_phenotype": ph,
                        "biobank_tier": t_name,
                        "n_blocks_in_tier": len(t_pts),
                        "mean_biobank_max_chisq": round(float(np.mean(t_bio_c)), 3),
                        "mean_lcl_max_chisq": round(float(np.mean(t_lcl_c)), 3),
                        "median_lcl_max_chisq": round(float(np.median(t_lcl_c)), 3),
                        "lcl_chisq_ratio_vs_bg": round(float(np.mean(t_lcl_c)) / max(1e-6, bg_lcl_mean), 4),
                        "frac_lcl_block_p_lt_1e-3": round(float(np.mean(t_lcl_p < 1e-3)), 4),
                        "tier_vs_bg_mannwhitney_p": round(float(p_mw), 5),
                        "global_spearman_rho_max_chisq": round(float(rho_max), 4),
                        "global_spearman_p_max_chisq": round(float(p_rho_max), 5),
                        "global_partial_spearman_rho_adj_nsnps": round(float(rho_max_adj), 4),
                        "global_partial_spearman_p_adj_nsnps": round(float(p_rho_max_adj), 5),
                        "global_spearman_rho_mean_chisq": round(float(rho_mean), 4),
                        "global_spearman_p_mean_chisq": round(float(p_rho_mean), 5),
                        "median_quantized_concordance_frac": round(concordant_above_below / n_tot, 4),
                        "median_quantized_binom_p": round(p_binom_med, 5),
                        "top10pct_overlap_obs_vs_exp": f"{obs_top10}/{round(exp_top10, 1)}",
                        "top10pct_fisher_or": round(float(or_top10), 3),
                        "top10pct_fisher_p": round(float(p_top10), 5),
                    })

    with open(f"{OUT_DIR}/biobank_unsigned_block_concordance.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(unsigned_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(unsigned_rows)

    # =========================================================================
    # ANALYSIS 2A & 2B: SIGNED QUANTIZED BLOCK TESTS (Both Options!)
    # =========================================================================
    print("Step 4: Running Analysis 2A (Allelic-Sign Concordance) & Analysis 2B (Germline TWAS Sign Concordance)...")
    # Load observed MAGE host-viral mRNA correlations (host_viral_mediation_20154.tsv)
    # Within-pop/batch correlation between BGLF5 and LMP1 in N=717 MAGE donors is r_bl = +0.5567
    R_BGLF5_LMP1 = 0.5567
    mage_rna_corr = {}  # geneSymbol -> {r_cn, r_bglf5, r_lmp1, q_bglf5, q_lmp1, r_lmp1_cond}
    for cand_p in [
        "/work/nfp8/EBV/metadata/host_viral_mediation_20154.tsv",
        f"{GWAS_DIR}/host_viral_mediation_20154.tsv",
        "/hpc/group/parrishlab/EBV_1kGP_MAGE/replication_2026/data/host_viral_mediation_20154.tsv",
    ]:
        if os.path.exists(cand_p):
            with open(cand_p) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    rb = float(row["r_bglf5"])
                    rl = float(row["r_lmp1"])
                    denom = math.sqrt(max(1e-12, (1.0 - rb * rb) * (1.0 - R_BGLF5_LMP1 * R_BGLF5_LMP1)))
                    rl_cond = (rl - rb * R_BGLF5_LMP1) / denom
                    mage_rna_corr[row["gene_symbol"]] = {
                        "r_cn": float(row["r_wgs_cn"]),
                        "r_bglf5": rb,
                        "q_bglf5": float(row["q_bglf5"]),
                        "r_lmp1": rl,
                        "q_lmp1": float(row["q_lmp1"]),
                        "r_lmp1_cond": rl_cond,
                        "q_lmp1_cond": float(row["q_lmp1"]),
                    }
            break
    print(f"  Loaded empirical MAGE host-viral mRNA correlations for {len(mage_rna_corr)} genes.")

    signed_rows = []

    for bio_name, bio_matched in [("Blood_UKB_AoU", blood_matched), ("Saliva_AoU", saliva_matched)]:
        for mask_mode in ["Unmasked_All", "Immune_MHC_Masked_SubThreshold"]:
            for ph in ["WGS_CN_ViralCond", "WGS_CN_PopAdj", "RNA_IELatRatio_ViralCond", "RNA_ProdFork_ViralCond", "RNA_BGLF5_ViralCond", "RNA_LMP1_ViralCond"]:

                # Deduplicate to 1 lead eQTL per unique variant (vid) and 1 lead eQTL per unique gene (highest PIP)
                block_best_eqtl = {}  # (chrom, b_idx) -> rec
                uniq_vid_eqtl = {}    # vid -> rec (highest PIP gene for that variant)
                uniq_gene_eqtl = {}   # gene -> rec (highest PIP variant for that gene)
                for g_row in eqtl_gene_rows:
                    vid = g_row["vid"]
                    if vid not in bio_matched or ph not in eqtl_lcl.get(vid, {}):
                        continue
                    chrom = g_row["chrom"]
                    b_idx = g_row["pos"] // 500000
                    is_masked = is_mhc_or_gwsig_block(chrom, b_idx)
                    if mask_mode == "Immune_MHC_Masked_SubThreshold" and is_masked:
                        continue
                    b_bio, se_bio, z_bio, p_bio, _ = bio_matched[vid]
                    z_lcl, p_lcl, b_lcl, se_lcl = eqtl_lcl[vid][ph]
                    afc = g_row["log2_afc"]
                    gene = g_row["gene"]
                    pip = g_row["pip"]
                    rec = (abs(z_bio), z_bio, p_bio, z_lcl, p_lcl, afc, gene, pip, is_masked)
                    if vid not in uniq_vid_eqtl or pip > uniq_vid_eqtl[vid][7]:
                        uniq_vid_eqtl[vid] = rec
                    if gene not in uniq_gene_eqtl or pip > uniq_gene_eqtl[gene][7]:
                        uniq_gene_eqtl[gene] = rec
                    cur = block_best_eqtl.get((chrom, b_idx))
                    if cur is None or abs(z_bio) > cur[0]:
                        block_best_eqtl[(chrom, b_idx)] = rec

                # --- Mode 2A_LCL_Manhattan_BlockLead: 5,457 LCL 500-kb block leads matched in Biobank ---
                lcl_lead_pts = []
                for (chrom, b_idx), ph_dict in lcl_blocks.items():
                    if ph not in ph_dict:
                        continue
                    is_masked = is_mhc_or_gwsig_block(chrom, b_idx)
                    if mask_mode == "Immune_MHC_Masked_SubThreshold" and is_masked:
                        continue
                    vid, pos, ref, alt, z_lcl, p_lcl, b_lcl, freq = ph_dict[ph]
                    if vid not in bio_matched:
                        continue
                    b_bio, se_bio, z_bio, p_bio, _ = bio_matched[vid]
                    lcl_lead_pts.append((abs(z_bio), z_bio, p_bio, z_lcl, p_lcl, 0.0, vid, 1.0, is_masked))

                for test_type, dataset in [
                    ("Option2A_BiobankLead_per_500kb_Block_AllelicSign", list(block_best_eqtl.values())),
                    ("Option2A_LCL_ManhattanLead_per_500kb_Block_AllelicSign", lcl_lead_pts),
                    ("Option2B_Unique_eQTL_CrossStudy_Sign_Unweighted", list(uniq_vid_eqtl.values())),
                ]:
                    tiers = defaultdict(list)
                    tiers["All_Included_Blocks_or_eQTLs"] = dataset
                    for item in dataset:
                        abs_z, z_bio, p_bio, z_lcl, p_lcl, afc, gene, pip, is_m = item
                        if p_bio < 5e-8:
                            tiers["Tier0_GWSig_P_lt_5e-8"].append(item)
                        elif p_bio < 1e-3:
                            tiers["Tier1_SubThreshold_P_lt_1e-3"].append(item)
                        elif p_bio < 0.05:
                            tiers["Tier2_Nominal_1e-3_to_0.05"].append(item)
                        else:
                            tiers["Tier3_Polygenic_Tail_P_ge_0.05"].append(item)
                        if p_bio < 0.05 and p_lcl < 0.05:
                            tiers["Tier_DualNominal_BioP05_and_LCLP05"].append(item)

                    for t_name, pts in sorted(tiers.items()):
                        if len(pts) < 5:
                            continue
                        z_b = np.array([x[1] for x in pts], dtype=np.float64)
                        z_l = np.array([x[3] for x in pts], dtype=np.float64)
                        # NOTE: Do NOT multiply both z_b and z_l by w_afc!
                        # Polarizing by sign(w_afc) on both sides preserves sign(z_b)*sign(z_l) without introducing w_afc^2 outlier leverage.
                        if test_type.startswith("Option2B"):
                            sgn_afc = np.array([1.0 if x[5] >= 0 else -1.0 for x in pts], dtype=np.float64)
                            z_b = sgn_afc * z_b
                            z_l = sgn_afc * z_l
                        non_zero = (z_b != 0) & (z_l != 0)
                        z_b = z_b[non_zero]
                        z_l = z_l[non_zero]
                        n_tot = len(z_b)
                        if n_tot < 5:
                            continue
                        n_same = int(np.sum(np.sign(z_b) == np.sign(z_l)))
                        frac_same = n_same / n_tot
                        p_binom = float(stats.binomtest(n_same, n_tot, 0.5).pvalue)
                        r_z, p_rz = stats.pearsonr(z_b, z_l)
                        rho_z, p_rhoz = stats.spearmanr(z_b, z_l)
                        slope, intercept, _, p_slope, se_slope = stats.linregress(z_b, z_l)

                        signed_rows.append({
                            "test_option": test_type,
                            "biobank_cohort": bio_name,
                            "mask_mode": mask_mode,
                            "lcl_phenotype": ph,
                            "biobank_tier": t_name,
                            "n_blocks_or_genes": n_tot,
                            "n_same_sign": n_same,
                            "sign_concordance_frac": round(frac_same, 4),
                            "binom_sign_p": round(p_binom, 5),
                            "z_pearson_r": round(float(r_z), 4),
                            "z_pearson_p": round(float(p_rz), 5),
                            "z_spearman_rho": round(float(rho_z), 4),
                            "z_spearman_p": round(float(p_rhoz), 5),
                            "z_slope_lcl_on_bio": round(float(slope), 5),
                            "z_slope_se": round(float(se_slope), 5),
                        })

                # --- Option 2B-2: Germline TWAS vs. Empirical MAGE mRNA Correlation Sign Test (Deduplicated 1 Lead eQTL per Unique Gene) ---
                rna_target = "r_bglf5" if "BGLF5" in ph or "IELat" in ph or "ProdFork" in ph else ("r_lmp1_cond" if "LMP1" in ph else "r_cn")
                twas_vs_rna_pts = []
                for gene, item in uniq_gene_eqtl.items():
                    abs_z, z_bio, p_bio, z_lcl, p_lcl, afc, _, pip, is_m = item
                    if gene not in mage_rna_corr or afc == 0:
                        continue
                    r_obs = mage_rna_corr[gene][rna_target]
                    if abs(r_obs) < 0.08:
                        continue
                    z_twas_bio = (1.0 if afc > 0 else -1.0) * z_bio
                    z_twas_lcl = (1.0 if afc > 0 else -1.0) * z_lcl
                    twas_vs_rna_pts.append((p_bio, p_lcl, z_twas_bio, z_twas_lcl, r_obs, pip))

                if len(twas_vs_rna_pts) >= 10:
                    for sub_label, z_idx, p_idx in [
                        ("Option2B_Biobank_TWAS_vs_Observed_MAGE_mRNA_ZeroOverlap", 2, 0),
                        ("Option2B_LCL_GWAS_TWAS_vs_MAGE_mRNA_InSampleOverlapBiased", 3, 1),
                    ]:
                        if sub_label.startswith("Option2B_LCL") and bio_name != "Blood_UKB_AoU":
                            continue
                        for sub_tier, filt_fn in [
                            ("All_Coupled_Genes_absR_ge_0.08", lambda x: True),
                            ("Nominal_GWAS_P_lt_0.05_and_absR_ge_0.08", lambda x: x[p_idx] < 0.05),
                            ("HighPIP_ge_0.5_and_Nominal_GWAS_P_lt_0.05", lambda x: x[p_idx] < 0.05 and x[5] >= 0.5),
                        ]:
                            sub_pts = [x for x in twas_vs_rna_pts if filt_fn(x)]
                            if len(sub_pts) < 5:
                                continue
                            zt = np.array([x[z_idx] for x in sub_pts], dtype=np.float64)
                            ro = np.array([x[4] for x in sub_pts], dtype=np.float64)
                            n_tot = len(zt)
                            n_same = int(np.sum(np.sign(zt) == np.sign(ro)))
                            p_bin = float(stats.binomtest(n_same, n_tot, 0.5).pvalue)
                            r_c, p_rc = stats.pearsonr(zt, ro)
                            rho_c, p_rhoc = stats.spearmanr(zt, ro)
                            signed_rows.append({
                                "test_option": sub_label,
                                "biobank_cohort": bio_name if "Biobank" in sub_label else "LCL_1kGP_MAGE_InSampleOverlap",
                                "mask_mode": mask_mode,
                                "lcl_phenotype": ph,
                                "biobank_tier": sub_tier,
                                "n_blocks_or_genes": n_tot,
                                "n_same_sign": n_same,
                                "sign_concordance_frac": round(n_same / n_tot, 4),
                                "binom_sign_p": round(p_bin, 5),
                                "z_pearson_r": round(float(r_c), 4),
                                "z_pearson_p": round(float(p_rc), 5),
                                "z_spearman_rho": round(float(rho_c), 4),
                                "z_spearman_p": round(float(p_rhoc), 5),
                                "z_slope_lcl_on_bio": 0.0,
                                "z_slope_se": 0.0,
                            })

    with open(f"{OUT_DIR}/biobank_signed_block_and_twas_concordance.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(signed_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(signed_rows)

    # =========================================================================
    # ANALYSIS 3: PARTITIONING HERITABILITY ACROSS HOST & VIRAL GENE PROGRAMS
    # =========================================================================
    print("Step 5: Running Analysis 3 (Partitioning Signal Across Host & Viral Gene Expression Programs)...")
    lareau_all_148 = set()
    lareau_non_hla_129 = set()
    lareau_non_chr6_12 = set()
    if os.path.exists(f"{LAREAU_DIR}/celltype-pathway-mapping/data/ebv_hits.txt"):
        with open(f"{LAREAU_DIR}/celltype-pathway-mapping/data/ebv_hits.txt") as f:
            lareau_all_148 = {line.strip() for line in f if line.strip()}
    if os.path.exists(f"{LAREAU_DIR}/celltype-pathway-mapping/data/EBV_0.0018_Ch6inc_HLArem.txt"):
        with open(f"{LAREAU_DIR}/celltype-pathway-mapping/data/EBV_0.0018_Ch6inc_HLArem.txt") as f:
            lareau_non_hla_129 = {line.strip() for line in f if line.strip()}
    if os.path.exists(f"{LAREAU_DIR}/celltype-pathway-mapping/data/EBV_0.0018_Ch6rem.txt"):
        with open(f"{LAREAU_DIR}/celltype-pathway-mapping/data/EBV_0.0018_Ch6rem.txt") as f:
            lareau_non_chr6_12 = {line.strip() for line in f if line.strip()}

    mhc_ag_genes = {
        "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
        "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DQB2",
        "HLA-DPA1", "HLA-DPB1", "TAP1", "TAP2", "TAPBP", "B2M", "CIITA", "NLRC5", "CD74",
        "ERAP1", "ERAP2", "LNPEP"
    }
    tnk_extrinsic_genes = lareau_non_chr6_12 | {
        "CTLA4", "CD27", "CD70", "SH2B3", "PTPN22", "ICOS", "CD226", "BACH2", "IFNG",
        "ITGAL", "LY9", "GSDMB", "NBEAL2", "SLC22A4", "FUT2", "EOMES", "IL2RA", "SLAMF1",
        "TNFRSF14", "CXCR5", "PRF1", "GZMB", "GZMA", "KLRK1", "KLRC1", "NKG7", "CD8A", "CD8B", "CD4"
    } - {"ERAP1", "ERAP2", "LNPEP"}

    epigenetic_silencer_genes = {
        "DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "UHRF1",
        "SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A",
        "TASOR", "MPHOSPH8", "PPHLN1", "MORC2", "SETDB1", "TRIM28", "CBX5",
        "PML", "SP100", "DAXX", "ATRX", "CTCF", "PARP1", "ZNF814"
    }
    bcell_lytic_tf_genes = {
        "CR2", "CD19", "CD81", "CD40", "RBPJ", "PAX5", "EBF1", "SPI1", "POU2F2", "POU2AF1",
        "IRF4", "IRF8", "NFKB1", "NFKB2", "REL", "RELA", "RELB", "TRAF1", "STAT3", "MYC",
        "PRDM1", "XBP1", "MEF2D", "KLF4", "SP1", "CREB1", "ATF1", "JUN", "FOS", "HIF1A", "MED21"
    }
    lmp1_sorted = sorted(mage_rna_corr.items(), key=lambda kv: abs(kv[1]["r_lmp1_cond"]), reverse=True)
    bglf5_sorted = sorted(mage_rna_corr.items(), key=lambda kv: abs(kv[1]["r_bglf5"]), reverse=True)
    mage_lmp1_top500 = {g for g, _ in lmp1_sorted[:500]}
    mage_bglf5_top500 = {g for g, _ in bglf5_sorted[:500]}

    def classify_gene_programs(gene, chrom, pos):
        progs = []
        is_mhc_region = (chrom == "6" and 25000000 <= pos <= 35000000)
        if gene in mhc_ag_genes or is_mhc_region or gene.startswith("HLA-"):
            progs.append("1_MHC_and_Antigen_Processing_HLA_ERAP")
        if gene in tnk_extrinsic_genes or (gene in lareau_all_148 and not is_mhc_region):
            progs.append("2_Biobank_T_NK_Cell_Extrinsic_Immune_Loci")
        if gene in epigenetic_silencer_genes or (gene.startswith("ZNF") and not is_mhc_region):
            progs.append("3_Epigenetic_CpG_Methylation_KZNF_Episome_Silencers")
        if gene in bcell_lytic_tf_genes:
            progs.append("4_BCell_Identity_and_Zp_Rp_Lytic_Switch_TFs")
        if gene in mage_lmp1_top500 and not is_mhc_region:
            progs.append("5_MAGE_LMP1_NFkB_Activated_Host_Program_Top500")
        if gene in mage_bglf5_top500 and not is_mhc_region:
            progs.append("6_MAGE_BGLF5_Lytic_Shutoff_Coupled_Program_Top500")
        if not progs and not is_mhc_region:
            progs.append("7_Background_NonMHC_Expressed_Genes")
        return progs

    # Deduplicate eqtl_gene_rows to 1 lead cis-eQTL per unique gene (highest PIP)
    lead_eqtl_per_gene = {}
    for g_row in eqtl_gene_rows:
        g = g_row["gene"]
        if g not in lead_eqtl_per_gene or g_row["pip"] > lead_eqtl_per_gene[g]["pip"]:
            lead_eqtl_per_gene[g] = g_row

    prog_buckets = defaultdict(lambda: defaultdict(list))  # prog -> trait_name -> list of (chi2, pval, b_key)
    for gene, g_row in lead_eqtl_per_gene.items():
        vid = g_row["vid"]
        chrom = g_row["chrom"]
        pos = g_row["pos"]
        b_key = (chrom, pos // 500000)
        progs = classify_gene_programs(gene, chrom, pos)
        if vid in blood_matched:
            _, _, zb, pb, _ = blood_matched[vid]
            for pr in progs:
                prog_buckets[pr]["Biobank_Blood_UKB_AoU"].append((zb * zb, pb, b_key))
        if vid in saliva_matched:
            _, _, zs, ps, _ = saliva_matched[vid]
            for pr in progs:
                prog_buckets[pr]["Biobank_Saliva_AoU"].append((zs * zs, ps, b_key))
        for ph in ["WGS_CN_ViralCond", "RNA_IELatRatio_ViralCond", "RNA_ProdFork_ViralCond", "RNA_BGLF5_ViralCond", "RNA_LMP1_ViralCond"]:
            if ph in eqtl_lcl.get(vid, {}):
                tl, pl, _, _ = eqtl_lcl[vid][ph]
                for pr in progs:
                    prog_buckets[pr][f"LCL_{ph}"].append((tl * tl, pl, b_key))

    traits_order = [
        "Biobank_Blood_UKB_AoU",
        "Biobank_Saliva_AoU",
        "LCL_WGS_CN_ViralCond",
        "LCL_RNA_IELatRatio_ViralCond",
        "LCL_RNA_ProdFork_ViralCond",
        "LCL_RNA_BGLF5_ViralCond",
        "LCL_RNA_LMP1_ViralCond",
    ]
    bg_means = {}
    bg_block_arrays = {}
    for tr in traits_order:
        bg_vals = [x[0] for x in prog_buckets["7_Background_NonMHC_Expressed_Genes"][tr]]
        bg_means[tr] = float(np.mean(bg_vals)) if bg_vals else 1.0
        b_map = defaultdict(list)
        for chi2_v, _, bk in prog_buckets["7_Background_NonMHC_Expressed_Genes"][tr]:
            b_map[bk].append(chi2_v)
        bg_block_arrays[tr] = np.array([np.mean(v) for v in b_map.values()], dtype=np.float64)

    prog_rows = []
    for pr in sorted(prog_buckets.keys()):
        for tr in traits_order:
            vals = prog_buckets[pr][tr]
            if not vals:
                continue
            c_arr = np.array([x[0] for x in vals], dtype=np.float64)
            p_arr = np.array([x[1] for x in vals], dtype=np.float64)
            bg_arr = np.array([x[0] for x in prog_buckets["7_Background_NonMHC_Expressed_Genes"][tr]], dtype=np.float64)
            _, p_mw_gene = stats.mannwhitneyu(c_arr, bg_arr, alternative="greater") if pr != "7_Background_NonMHC_Expressed_Genes" else (0.0, 1.0)
            # Collapse tandem gene clusters by 500-kb LD block
            pr_bmap = defaultdict(list)
            for chi2_v, _, bk in vals:
                pr_bmap[bk].append(chi2_v)
            pr_blk_arr = np.array([np.mean(v) for v in pr_bmap.values()], dtype=np.float64)
            _, p_mw_blk = stats.mannwhitneyu(pr_blk_arr, bg_block_arrays[tr], alternative="greater") if pr != "7_Background_NonMHC_Expressed_Genes" and len(pr_blk_arr) >= 3 else (0.0, 1.0)
            prog_rows.append({
                "gene_program": pr,
                "trait": tr,
                "n_unique_genes": len(vals),
                "n_500kb_ld_blocks": len(pr_blk_arr),
                "mean_chisq": round(float(np.mean(c_arr)), 3),
                "median_chisq": round(float(np.median(c_arr)), 3),
                "max_chisq": round(float(np.max(c_arr)), 2),
                "min_p": f"{float(np.min(p_arr)):.3e}",
                "fold_enrichment_vs_bg": round(float(np.mean(c_arr)) / max(1e-6, bg_means[tr]), 3),
                "frac_p_lt_0.01": round(float(np.mean(p_arr < 0.01)), 4),
                "frac_p_lt_0.001": round(float(np.mean(p_arr < 0.001)), 4),
                "gene_mannwhitney_enrichment_p": f"{float(p_mw_gene):.3e}",
                "block_mannwhitney_enrichment_p": f"{float(p_mw_blk):.3e}",
            })

    with open(f"{OUT_DIR}/biobank_lcl_gene_program_partition.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prog_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(prog_rows)

    # =========================================================================
    # ANALYSIS 4A: BIOBANK HIGH-PRIOR LOOKUP TABLE (Exact Lead SNP + +/-100kb Window)
    # =========================================================================
    print("Step 6: Running Analysis 4A (Biobank High-Prior Lookup Table via remote tabix + LCL windows)...")
    target_snps = []
    for g, ch, st, en, cat in PRIORITY_LOCI:
        if g in blood_loci:
            pval, chi2, chrom, pos, a1, a2, beta, se, z = blood_loci[g]
            target_snps.append((chrom, pos, a1, a2))
        if g in saliva_loci:
            pval, chi2, chrom, pos, a1, a2, beta, se, z = saliva_loci[g]
            target_snps.append((chrom, pos, a1, a2))

    exact_lcl_at_bio_leads = test_exact_snps_in_lcls(target_snps)
    print(f"  Tested {len(exact_lcl_at_bio_leads)} exact Biobank locus lead SNPs in 1kGP LCL genotypes via tabix (FWL covariate-residualized).")

    cand55_lcl = defaultdict(dict)  # gene -> {ph: (min_p, best_vid, best_t)}
    suggestive_by_chr = defaultdict(lambda: defaultdict(list))
    for c in range(1, 23):
        chrom = str(c)
        for ph in LCL_PHENOS:
            fpath = f"{CHR_DIR}/chr{chrom}_{ph}_candidate55.tsv"
            if os.path.exists(fpath):
                with open(fpath) as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        cand55_lcl[row["gene"]][ph] = (float(row["p"]), row["id"], float(row["t_stat"]))
            spath = f"{CHR_DIR}/chr{chrom}_{ph}_suggestive.tsv"
            if os.path.exists(spath):
                with open(spath) as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        suggestive_by_chr[chrom][ph].append((int(row["pos"]), float(row["p"]), row["id"], float(row["t_stat"])))

    for g, ch, st, en, cat in PRIORITY_LOCI:
        w_st, w_en = max(1, st - 100000), en + 100000
        for ph in LCL_PHENOS:
            best_p, best_id, best_t = cand55_lcl[g].get(ph, (1.0, "NA", 0.0))
            for pos, pl, vid, tl in suggestive_by_chr[ch].get(ph, []):
                if w_st <= pos <= w_en and pl < best_p:
                    best_p, best_id, best_t = pl, vid, tl
            for g_row in eqtl_gene_rows:
                if g_row["gene"] == g or (g_row["chrom"] == ch and w_st <= g_row["pos"] <= w_en):
                    vid = g_row["vid"]
                    if ph in eqtl_lcl.get(vid, {}):
                        tl, pl, _, _ = eqtl_lcl[vid][ph]
                        if pl < best_p:
                            best_p, best_id, best_t = pl, vid, tl
            for b_idx in range(w_st // 500000, (w_en // 500000) + 1):
                if ph in lcl_blocks.get((ch, b_idx), {}):
                    vid, pos, _, _, tl, pl, _, _ = lcl_blocks[(ch, b_idx)][ph]
                    if w_st <= pos <= w_en and pl < best_p:
                        best_p, best_id, best_t = pl, vid, tl
            cand55_lcl[g][ph] = (best_p, best_id, best_t)

    lookup_rows = []
    for g, ch, st, en, cat in PRIORITY_LOCI:
        b_rec = blood_loci.get(g)
        s_rec = saliva_loci.get(g)
        b_p = b_rec[0] if b_rec else 1.0
        b_z = b_rec[8] if b_rec else 0.0
        b_snp = f"chr{b_rec[2]}:{b_rec[3]}:{b_rec[4]}:{b_rec[5]}" if b_rec else "NA"
        s_p = s_rec[0] if s_rec else 1.0
        s_z = s_rec[8] if s_rec else 0.0
        s_snp = f"chr{s_rec[2]}:{s_rec[3]}:{s_rec[4]}:{s_rec[5]}" if s_rec else "NA"

        # Exact LCL stats at the Biobank Blood lead SNP ("NA" if variant not present/matched in 1kGP)
        ex_cn_z_str, ex_cn_p_str = "NA", "NA"
        ex_ie_z_str, ex_ie_p_str = "NA", "NA"
        if b_rec and (b_rec[2], b_rec[3]) in exact_lcl_at_bio_leads:
            ex_dict = exact_lcl_at_bio_leads[(b_rec[2], b_rec[3])]
            if "WGS_CN_ViralCond" in ex_dict:
                ex_cn_z_str = f"{ex_dict['WGS_CN_ViralCond'][2]:.2f}"
                ex_cn_p_str = f"{ex_dict['WGS_CN_ViralCond'][3]:.3g}"
            if "RNA_IELatRatio_ViralCond" in ex_dict:
                ex_ie_z_str = f"{ex_dict['RNA_IELatRatio_ViralCond'][2]:.2f}"
                ex_ie_p_str = f"{ex_dict['RNA_IELatRatio_ViralCond'][3]:.3g}"

        rna_info = mage_rna_corr.get(g, {})
        lookup_rows.append({
            "gene": g,
            "chrom": ch,
            "window_grch38": f"chr{ch}:{st}-{en}",
            "biological_compartment": cat,
            "biobank_blood_lead_snp": b_snp,
            "biobank_blood_z": round(float(b_z), 2),
            "biobank_blood_p": f"{float(b_p):.2e}",
            "biobank_saliva_lead_snp": s_snp,
            "biobank_saliva_z": round(float(s_z), 2),
            "biobank_saliva_p": f"{float(s_p):.2e}",
            "lcl_exact_blood_lead_CN_z": ex_cn_z_str,
            "lcl_exact_blood_lead_CN_p": ex_cn_p_str,
            "lcl_exact_blood_lead_IELat_z": ex_ie_z_str,
            "lcl_exact_blood_lead_IELat_p": ex_ie_p_str,
            "lcl_window_minP_WGS_CN_ViralCond": f"{cand55_lcl[g]['WGS_CN_ViralCond'][0]:.2e}" if cand55_lcl[g]['WGS_CN_ViralCond'][1] != "NA" else "NA",
            "lcl_window_minP_RNA_IELatRatio_ViralCond": f"{cand55_lcl[g]['RNA_IELatRatio_ViralCond'][0]:.2e}" if cand55_lcl[g]['RNA_IELatRatio_ViralCond'][1] != "NA" else "NA",
            "lcl_window_minP_RNA_ProdFork_ViralCond": f"{cand55_lcl[g]['RNA_ProdFork_ViralCond'][0]:.2e}" if cand55_lcl[g]['RNA_ProdFork_ViralCond'][1] != "NA" else "NA",
            "lcl_window_minP_RNA_LMP1_ViralCond": f"{cand55_lcl[g]['RNA_LMP1_ViralCond'][0]:.2e}" if cand55_lcl[g]['RNA_LMP1_ViralCond'][1] != "NA" else "NA",
            "mage_host_mRNA_r_BGLF5": round(float(rna_info.get("r_bglf5", 0.0)), 3) if rna_info else "NA",
            "mage_host_mRNA_r_LMP1_cond": round(float(rna_info.get("r_lmp1_cond", 0.0)), 3) if rna_info else "NA",
        })

    with open(f"{OUT_DIR}/biobank_ebv_loci_comparison.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lookup_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(lookup_rows)

    # =========================================================================
    # ANALYSIS 4B: IN VIVO BIOBANK (UKB + AoU) vs. 1kGP LCL VIRAL VARIANT COMPARISON
    # =========================================================================
    print("Step 7: Running Analysis 4B (Comparing 1kGP LCL Viral Variants vs. UKB & AoU In Vivo Blood Pileups)...")
    ukb_ebv = {}
    aou_ebv = {}
    for tag, fpath, store in [
        ("UKB", f"{LAREAU_DIR}/viral-sequences/data/UKB_chrEBV.tsv.gz", ukb_ebv),
        ("AoU", f"{LAREAU_DIR}/viral-sequences/data/AOU_chrEBV.tsv.gz", aou_ebv),
    ]:
        if os.path.exists(fpath):
            with gzip.open(fpath, "rt") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    bp = int(row["BP"])
                    tot = int(row["TOTAL"])
                    store[bp] = {
                        "REF": row["REF"],
                        "TOTAL": tot,
                        "A": int(row["A"]),
                        "C": int(row["C"]),
                        "G": int(row["G"]),
                        "T": int(row["T"]),
                    }

    # Load our 112 1kGP viral variants from viral_qtl_scan_112.tsv
    viral_rows = []
    for cand_vp in [
        "/work/nfp8/EBV/metadata/viral_qtl_scan_112.tsv",
        f"{GWAS_DIR}/viral_qtl_scan_112.tsv",
        "/hpc/group/parrishlab/EBV_1kGP_MAGE/replication_2026/data/viral_qtl_scan_112.tsv",
    ]:
        if os.path.exists(cand_vp):
            with open(cand_vp) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    pos = int(row["pos"])
                    ref = row["ref"]
                    alt = row["alt"]
                    n_wgs = int(row["n_wgs_2449"])
                    lcl_af = n_wgs / 2449.0
                    r_cn = float(row["r_cn_within_pop_2449"])
                    t_cn = r_cn * math.sqrt(2447.0 / max(1e-12, 1.0 - r_cn * r_cn))
                    p_cn = float(2.0 * stats.t.sf(abs(t_cn), df=2447))
                    u_rec = ukb_ebv.get(pos, {})
                    a_rec = aou_ebv.get(pos, {})
                    u_tot = u_rec.get("TOTAL", 0)
                    a_tot = a_rec.get("TOTAL", 0)
                    is_true_snv = (len(ref) == 1 and len(alt) == 1 and alt in ("A", "C", "G", "T"))
                    if is_true_snv:
                        u_alt_ct = u_rec.get(alt, 0)
                        a_alt_ct = a_rec.get(alt, 0)
                        u_af = u_alt_ct / u_tot if u_tot > 0 else 0.0
                        a_af = a_alt_ct / a_tot if a_tot > 0 else 0.0
                        mean_invivo_af = 0.5 * (u_af + a_af)
                        origin_class = (
                            "WildType_Polymorphic_InVivo_and_LCL_SNV"
                            if mean_invivo_af >= 0.05
                            else ("Rare_InVivo_0.5_to_5pct_SNV" if mean_invivo_af >= 0.005 else "B95-8_Lab_Culture_Enriched_SNV_lt_0.5pct_InVivo")
                        )
                        u_alt_str = str(u_alt_ct)
                        u_af_str = f"{u_af:.4f}"
                        a_alt_str = str(a_alt_ct)
                        a_af_str = f"{a_af:.4f}"
                    else:
                        # Single-base A/C/G/T pileups in UKB_chrEBV.tsv.gz do NOT assay indels
                        origin_class = "Indel_NotAssayedInSingleBasePileup"
                        u_alt_str = "NA_Indel"
                        u_af_str = "NA_Indel"
                        a_alt_str = "NA_Indel"
                        a_af_str = "NA_Indel"

                    viral_rows.append({
                        "variant_id": row["site"],
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "variant_type": "SNV" if is_true_snv else "INDEL",
                        "viral_gene": row.get("gene", ""),
                        "kinetic_class": row.get("kinetic_class", ""),
                        "lcl_1kgp_carrier_freq": round(lcl_af, 4),
                        "lcl_wgs_cn_r_within_pop": round(r_cn, 4),
                        "lcl_wgs_cn_p": f"{p_cn:.3e}",
                        "lcl_ie_rna_r": round(float(row["r_ie_within_pop_batch_717"]), 4),
                        "ukb_invivo_total_reads": u_tot,
                        "ukb_invivo_alt_reads": u_alt_str,
                        "ukb_invivo_alt_freq": u_af_str,
                        "aou_invivo_total_reads": a_tot,
                        "aou_invivo_alt_reads": a_alt_str,
                        "aou_invivo_alt_freq": a_af_str,
                        "origin_classification": origin_class,
                    })
            break

    if viral_rows:
        with open(f"{OUT_DIR}/biobank_vs_1kgp_viral_variants.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(viral_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(viral_rows)

    print("All 4 Biobank vs. LCL cross-study analyses completed successfully!")


if __name__ == "__main__":
    main()
