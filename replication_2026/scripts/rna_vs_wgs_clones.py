#!/usr/bin/env python3
"""Compare WGS BCR clonotypes (PoN-filtered: valid V+J and private in <=2 donors across 2,502)
against MAGE RNA-seq TRUST4 BCR clonotypes without any artificial count<=20 cap."""
import csv, glob, os, statistics as st
from collections import Counter, defaultdict

BASE = '/work/nfp8/EBV'
wgs_raw = {}
cdr3_donors = defaultdict(set)
for f in glob.glob(f'{BASE}/bcr/wgs/*_report.tsv'):
    s = os.path.basename(f).replace('_report.tsv', '')
    clones = []
    for line in open(f):
        if line.startswith('#') or line.startswith('count'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 8 or not p[4].startswith('IGH'): continue
        cnt = int(p[0]); cdr3 = p[3]; v = p[4]; j = p[6]
        if v == '*' or j == '*' or len(cdr3) < 8 or not cdr3.startswith('C'): continue
        clones.append((cnt, cdr3, v, j))
        cdr3_donors[cdr3].add(s)
    wgs_raw[s] = clones

public_pon = {c for c, ds in cdr3_donors.items() if len(ds) > 2}

rna_map = defaultdict(list)
for r in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sequencing_library.metadata.MAGE.v1.0.txt'), delimiter='\t'):
    run = r['SRA_accession']; s = r['sample_kgpID']
    rep = f'{BASE}/trust4/{run}/{run}_report.tsv'
    if os.path.exists(rep): rna_map[s].append(rep)

rna_clones = {}
for s, reps in rna_map.items():
    c_cnt = Counter()
    for rep in reps:
        for line in open(rep):
            if line.startswith('#') or line.startswith('count'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) >= 8 and p[4].startswith('IGH') and len(p[3]) >= 8:
                c_cnt[p[3]] += int(p[0])
    if c_cnt: rna_clones[s] = c_cnt

paired = [s for s in rna_clones if s in wgs_raw]
rec_rates = []
top1_in_rna = top1_in_top5 = top1_is_top1 = 0
n_eval = 0
for s in paired:
    wc = sorted([(cnt, cdr3) for cnt, cdr3, v, j in wgs_raw[s] if cdr3 not in public_pon], reverse=True)
    if not wc: continue
    n_eval += 1
    rc = rna_clones[s]
    r_top5 = {c for c, _ in rc.most_common(5)}
    r_top1 = rc.most_common(1)[0][0]
    rec_rates.append(sum(1 for _, c in wc if c in rc) / len(wc))
    t1 = wc[0][1]
    if t1 in rc: top1_in_rna += 1
    if t1 in r_top5: top1_in_top5 += 1
    if t1 == r_top1: top1_is_top1 += 1

print(f'Paired donors (PoN-cleaned, no count<=20 cap): n = {n_eval}')
print(f'  Mean WGS clones recovered in RNA : {st.mean(rec_rates)*100:.1f}% (median {st.median(rec_rates)*100:.1f}%)')
print(f'  Top WGS clone present in RNA     : {top1_in_rna}/{n_eval} ({top1_in_rna/n_eval*100:.1f}%)')
print(f'  Top WGS clone in RNA top 5       : {top1_in_top5}/{n_eval} ({top1_in_top5/n_eval*100:.1f}%)')
print(f'  Top WGS clone is top RNA clone   : {top1_is_top1}/{n_eval} ({top1_is_top1/n_eval*100:.1f}%)')
