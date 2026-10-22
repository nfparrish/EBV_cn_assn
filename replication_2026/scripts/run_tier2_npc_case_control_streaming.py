#!/usr/bin/env python3
"""
run_tier2_npc_case_control_streaming.py
=======================================
Zero-disk HTTP streaming and within-study case-control association pipeline
for testing EBV Copy-Number QTLs in primary human clinical cohorts containing
internal healthy controls:

  1. PRJNA1268153 (SYSUCC Two-Cohort NPC Case-Control EBV WGS, N = 732 runs):
       - 202 NPC Cases (_ca_) vs. 530 Healthy Controls (_co_)
       - Internal discovery (cohort1) and replication (cohort2) split
  2. PRJNA522388 (Xu et al. 2019 Nature Genetics, SYSUCC, N = 270 runs):
       - 50 Endemic Healthy Saliva Controls (HS) (+ 4 non-endemic NHS)
       - 53 NPC Patient Saliva (NPCS) -> compartment-matched Saliva vs. Saliva
       - 113 Endemic NPC Tumor Biopsies (NPCT) -> Tumor vs. Healthy Saliva

Strictly enforces WITHIN-STUDY association testing (no cross-study pooling)
to prevent geographic EBV lineage and hybrid-capture batch confounding.

Usage:
  python3 run_tier2_npc_case_control_streaming.py --study PRJNA1268153 --workers 12 --mb-per-run 8
  python3 run_tier2_npc_case_control_streaming.py --study PRJNA522388  --workers 12 --full-stream
"""

import argparse
import concurrent.futures
import csv
import json
import math
import os
import ssl
import sys
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

import numpy as np
from scipy import stats
import statsmodels.api as sm

CTX = ssl._create_unverified_context()
RC_MAP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]


def fetch_ebv_reference() -> str:
    """Fetch NC_007605.1 (chrEBV, 171,823 bp) sequence from NCBI E-utilities."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=NC_007605.1&rettype=fasta&retmode=text"
    )
    req = urllib.request.urlopen(url, timeout=30, context=CTX)
    fa = req.read().decode("utf-8")
    return "".join(l.strip() for l in fa.splitlines() if not l.startswith(">")).upper()


def build_probes(ref: str):
    """
    Build exact 25-mer forward/reverse-complement probes for:
      - 12,129:T>A (IR1 / W2-W1 intron / #1 LCL CN QTL)
      - 12,024:G>T (IR1 Copy 0 boundary tag)
      - 36,480:T>C (EBNA-2 p.Ser89Pro in N-terminal Poly-Proline Repeat)
      - 36,494:G>C (EBNA-2 p.Pro93= in N-terminal Poly-Proline Repeat)
      - 162,476:C>T & 163,364:C>T (BALF2 East Asian high-risk NPC lineage tags)
      - 20 probes across the 3,072 bp IR1 (BamHI-W) repeat unit (12,200..14,900)
      - 40 probes across single-copy unique EBV regions (45,000..160,000)
    """
    snp_defs = [
        ("12129_T_A", 12129, "T", "A"),
        ("12024_G_T", 12024, "G", "T"),
        ("36480_T_C", 36480, "T", "C"),
        ("36494_G_C", 36494, "G", "C"),
        ("162476_C_T", 162476, "C", "T"),
        ("163364_C_T", 163364, "C", "T"),
    ]
    snp_probes = {}
    for name, pos, r_base, a_base in snp_defs:
        idx = pos - 1
        left = ref[idx - 12 : idx]
        right = ref[idx + 1 : idx + 13]
        rk = left + r_base + right
        ak = left + a_base + right
        snp_probes[name] = (rk, revcomp(rk), ak, revcomp(ak))

    ir1_probes = []
    for p in range(12200, 14900, 135):
        k = ref[p : p + 25]
        ir1_probes.append((k, revcomp(k)))

    uniq_probes = []
    for p in range(45000, 160000, 2875):
        k = ref[p : p + 25]
        uniq_probes.append((k, revcomp(k)))

    return snp_probes, ir1_probes, uniq_probes


def fetch_study_runs(study_acc: str):
    """Fetch run metadata and FASTQ FTP links from ENA Portal API."""
    query = f'study_accession="{study_acc}" OR secondary_study_accession="{study_acc}"'
    fields = (
        "run_accession,sample_accession,sample_title,sample_alias,"
        "run_alias,isolation_source,fastq_bytes,fastq_ftp"
    )
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query="
        + urllib.parse.quote(query)
        + f"&fields={fields}&format=json&limit=2000"
    )
    req = urllib.request.urlopen(url, timeout=30, context=CTX)
    return json.loads(req.read().decode("utf-8"), strict=False)


def classify_sample(study_acc: str, run: dict):
    """Assign strictly internal case/control and compartment labels within each study."""
    s_alias = run.get("sample_alias", "")
    r_alias = run.get("run_alias", "")
    s_title = run.get("sample_title", "")

    if study_acc == "PRJNA1268153":
        cohort = "cohort1" if "cohort1" in r_alias else ("cohort2" if "cohort2" in r_alias else "unknown")
        if "_ca_" in r_alias:
            return {"group": "NPC_Case", "is_case": 1, "subcohort": cohort, "compartment": "swab_or_saliva"}
        elif "_co_" in r_alias:
            return {"group": "Healthy_Control", "is_case": 0, "subcohort": cohort, "compartment": "swab_or_saliva"}
        return None

    if study_acc == "PRJNA522388":
        tag = s_alias.strip()
        if tag.startswith("HS"):
            return {"group": "Healthy_Saliva", "is_case": 0, "subcohort": "endemic", "compartment": "saliva"}
        elif tag.startswith("NPCS"):
            return {"group": "NPC_Saliva", "is_case": 1, "subcohort": "endemic", "compartment": "saliva"}
        elif tag.startswith("NPCT"):
            return {"group": "NPC_Tumor", "is_case": 1, "subcohort": "endemic", "compartment": "tumor"}
        return None

    return None


def stream_and_quantify_run(run: dict, snp_probes, ir1_probes, uniq_probes, max_bytes: int, full_stream: bool):
    """Stream FASTQ R1 over HTTPS in memory and count exact 25-mer matches."""
    ftp_field = run.get("fastq_ftp", "")
    if not ftp_field:
        return None
    r1_ftp = ftp_field.split(";")[0]
    url = "https://" + r1_ftp

    headers = {}
    if not full_stream and max_bytes > 0:
        headers["Range"] = f"bytes=0-{max_bytes - 1}"

    req = urllib.request.Request(url, headers=headers)
    decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
    snp_counts = {k: [0, 0] for k in snp_probes}
    ir1_counts = [0] * len(ir1_probes)
    uniq_counts = [0] * len(uniq_probes)
    n_reads = 0

    try:
        with urllib.request.urlopen(req, timeout=45, context=CTX) as resp:
            leftover = ""
            bytes_read = 0
            while True:
                chunk = resp.read(262144)
                if not chunk:
                    break
                bytes_read += len(chunk)
                try:
                    txt = leftover + decomp.decompress(chunk).decode("ascii", errors="ignore")
                except Exception:
                    break
                lines = txt.split("\n")
                leftover = lines.pop()
                # Process sequence lines (line index 1 mod 4)
                for i in range(1, len(lines), 4):
                    rd = lines[i]
                    n_reads += 1
                    for name, (rf, rr, af, ar) in snp_probes.items():
                        if rf in rd or rr in rd:
                            snp_counts[name][0] += 1
                        if af in rd or ar in rd:
                            snp_counts[name][1] += 1
                    for idx, (kf, kr) in enumerate(ir1_probes):
                        if kf in rd or kr in rd:
                            ir1_counts[idx] += 1
                    for idx, (kf, kr) in enumerate(uniq_probes):
                        if kf in rd or kr in rd:
                            uniq_counts[idx] += 1
                if not full_stream and max_bytes > 0 and bytes_read >= max_bytes:
                    break
    except Exception as exc:
        return {"run_accession": run.get("run_accession"), "error": str(exc)}

    mean_ir1 = float(np.mean(ir1_counts))
    mean_uniq = float(np.mean(uniq_counts))
    tot_uniq = int(np.sum(uniq_counts))
    ir1_cn_ratio = mean_ir1 / mean_uniq if mean_uniq > 0 else np.nan

    out = {
        "run_accession": run.get("run_accession"),
        "sample_alias": run.get("sample_alias") or run.get("run_alias"),
        "n_reads_streamed": n_reads,
        "tot_uniq_hits": tot_uniq,
        "mean_ir1_depth": round(mean_ir1, 3),
        "mean_uniq_depth": round(mean_uniq, 3),
        "ir1_cn_ratio": round(ir1_cn_ratio, 4) if not np.isnan(ir1_cn_ratio) else "NA",
    }
    for name, (rcnt, acnt) in snp_counts.items():
        tot = rcnt + acnt
        af = acnt / tot if tot > 0 else np.nan
        out[f"{name}_ref"] = rcnt
        out[f"{name}_alt"] = acnt
        out[f"{name}_af"] = round(af, 4) if not np.isnan(af) else "NA"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", required=True, choices=["PRJNA1268153", "PRJNA522388"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--mb-per-run", type=int, default=8)
    ap.add_argument("--full-stream", action="store_true")
    ap.add_argument("--out-tsv", type=str, default=None)
    args = ap.parse_args()

    ref = fetch_ebv_reference()
    snp_probes, ir1_probes, uniq_probes = build_probes(ref)
    runs = fetch_study_runs(args.study)

    labeled_runs = []
    for r in runs:
        meta = classify_sample(args.study, r)
        if meta is not None:
            r.update(meta)
            labeled_runs.append(r)

    print(f"[{args.study}] Found {len(labeled_runs)} case/control labeled runs.")
    max_bytes = args.mb_per_run * 1024 * 1024
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(
                stream_and_quantify_run, r, snp_probes, ir1_probes, uniq_probes, max_bytes, args.full_stream
            ): r
            for r in labeled_runs
        }
        for fut in concurrent.futures.as_completed(futs):
            r = futs[fut]
            res = fut.result()
            if res and "error" not in res:
                res.update({
                    "group": r["group"],
                    "is_case": r["is_case"],
                    "subcohort": r["subcohort"],
                    "compartment": r["compartment"],
                })
                results.append(res)

    out_path = args.out_tsv or f"tier2_{args.study}_streaming_results.tsv"
    if results:
        keys = list(results[0].keys())
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(results)
        print(f"Wrote {len(results)} rows to {out_path}")


if __name__ == "__main__":
    main()
