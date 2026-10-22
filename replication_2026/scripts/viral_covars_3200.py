#!/usr/bin/env python3
"""
Extract viral IR1 repeat copy number and key viral QTL / strain-barcode genotypes
across all 3,200 NYGC 30x 1kGP samples (2,502 unrelateds + 698 trio relatives)
directly from local `/work/nfp8/EBV/slices/{sample}.chrEBV.bam`.
"""
import os
import subprocess
import sys
import pysam

BASE = "/work/nfp8/EBV"
SAMTOOLS = "/hpc/group/parrishlab/envs/samtools-env/bin/samtools"
OUT_DIR = f"{BASE}/host_gwas/viral_covars_chunks"

# Key viral positions (1-based NC_007605.1 coordinates)
# Includes replicated within-pop CN-QTLs (12129:T>A, 163858:G>T) and top strain barcodes
VIRAL_SITES = [
    ("v12129_TA", 12129, "T", "A"),
    ("v163858_GT", 163858, "G", "T"),
    ("v91192_GA", 91192, "G", "A"),
    ("v132965_GGC", 132965, "G", "GC"),
    ("v162476_TC", 162476, "T", "C"),
    ("v169200_CT", 169200, "C", "T"),
    ("v38566_GA", 38566, "G", "A"),
    ("v65013_CT", 65013, "C", "T"),
    ("v81920_CT", 81920, "C", "T"),
    ("v103742_AG", 103742, "A", "G"),
    ("v147913_GA", 147913, "G", "A"),
]


def get_ir1_metrics(sample, bam_path):
    ir1_file = f"{BASE}/ir1/{sample}.tsv"
    if os.path.exists(ir1_file):
        with open(ir1_file) as f:
            lines = [l.strip().split("\t") for l in f if l.strip()]
            if len(lines) >= 2 and len(lines[1]) >= 6 and lines[1][4] != "NA":
                ir1_d = float(lines[1][1])
                cp_d = float(lines[1][2])
                uniq_d = float(lines[1][3])
                ratio_u = float(lines[1][4])
                return ir1_d, cp_d, uniq_d, ratio_u * (22737.0 / 3072.0)

    # Compute directly from chrEBV.bam using samtools depth -a -Q 0
    cmd = f"{SAMTOOLS} depth -a -Q 0 -r chrEBV {bam_path}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    i_sum = i_n = c_sum = c_n = u_sum = u_n = 0
    for line in proc.stdout:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        p = int(parts[1])
        d = int(parts[2])
        if 12301 <= p <= 35037:
            i_sum += d
            i_n += 1
        elif 10000 <= p <= 11500:
            c_sum += d
            c_n += 1
        if (p < 12301) or (35037 < p < 38382) or (39733 < p < 95984) or (96575 < p < 139724) or (151551 < p < 169740):
            u_sum += d
            u_n += 1
    proc.wait()
    ir1_d = i_sum / i_n if i_n else 0.0
    cp_d = c_sum / c_n if c_n else 0.0
    uniq_d = u_sum / u_n if u_n else 0.0
    ir1_copies = (ir1_d / uniq_d) * (22737.0 / 3072.0) if uniq_d > 0 else 0.0
    # Also cache in /work/nfp8/EBV/ir1/{sample}.tsv
    try:
        with open(ir1_file, "w") as out:
            out.write("sample\tir1_depth_mapq0\tcp_depth\tunique_depth\tir1_copies_vs_unique\tir1_copies_vs_Cp\n")
            ru = f"{ir1_d / uniq_d:.3f}" if uniq_d > 0 else "NA"
            rc = f"{ir1_d / cp_d:.3f}" if cp_d > 0 else "NA"
            out.write(f"{sample}\t{ir1_d:.4f}\t{cp_d:.4f}\t{uniq_d:.4f}\t{ru}\t{rc}\n")
    except Exception:
        pass
    return ir1_d, cp_d, uniq_d, ir1_copies


def get_site_vaf(bam, pos_1based, ref, alt):
    pos0 = pos_1based - 1
    is_ins = len(alt) > len(ref)
    tot = 0
    alt_cnt = 0
    for col in bam.pileup("chrEBV", pos0, pos0 + 1, truncate=True, min_mapping_quality=20, min_base_quality=20, max_depth=50000):
        if col.pos != pos0:
            continue
        for pr in col.pileups:
            if pr.is_refskip:
                continue
            tot += 1
            if is_ins:
                if pr.indel == len(alt) - len(ref):
                    alt_cnt += 1
            else:
                if not pr.is_del and pr.query_position is not None:
                    base = pr.alignment.query_sequence[pr.query_position].upper()
                    if base == alt:
                        alt_cnt += 1
    if tot < 5:
        return "NA"
    return f"{alt_cnt / tot:.4f}"


def main():
    task_id = int(sys.argv[1])
    num_tasks = int(sys.argv[2])
    os.makedirs(OUT_DIR, exist_ok=True)

    samples = []
    with open(f"{BASE}/metadata/all_cn_3200.tsv") as f:
        for line in f:
            p = line.strip().split("\t")
            if p:
                samples.append(p[0])
    samples = sorted(set(samples))
    chunk = [s for idx, s in enumerate(samples) if (idx % num_tasks) == (task_id - 1)]

    out_path = f"{OUT_DIR}/chunk_{task_id:02d}.tsv"
    header = ["sample", "ir1_depth_mapq0", "cp_depth", "unique_depth", "ir1_copies"] + [s[0] for s in VIRAL_SITES]
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for s in chunk:
            bam_path = f"{BASE}/slices/{s}.chrEBV.bam"
            if not os.path.exists(bam_path):
                continue
            ir1_d, cp_d, uniq_d, ir1_copies = get_ir1_metrics(s, bam_path)
            bam = pysam.AlignmentFile(bam_path, "rb")
            vafs = [get_site_vaf(bam, pos, ref, alt) for _, pos, ref, alt in VIRAL_SITES]
            bam.close()
            row = [s, f"{ir1_d:.2f}", f"{cp_d:.2f}", f"{uniq_d:.2f}", f"{ir1_copies:.4f}"] + vafs
            out.write("\t".join(row) + "\n")
    print(f"Task {task_id}/{num_tasks} finished ({len(chunk)} samples) -> {out_path}")


if __name__ == "__main__":
    main()
