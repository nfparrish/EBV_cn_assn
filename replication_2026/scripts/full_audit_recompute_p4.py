import csv, glob, math, os, random, statistics as st
from collections import defaultdict, Counter

BASE = '/work/nfp8/EBV'
RU = 22737 / 3072.0

pop, spop = {}, {}
for line in open(f'{BASE}/metadata/ped3202.txt'):
    f = line.split()
    if len(f) >= 7 and f[1] != 'SampleID':
        pop[f[1]] = f[5]; spop[f[1]] = f[6]

cn, wt = {}, {}
for line in open(f'{BASE}/metadata/all_cn_3200.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 10:
        try:
            v = float(f[3]); w = float(f[9])
            if v > 0: cn[f[0]] = v; wt[f[0]] = w
        except ValueError: pass

unrel_2502 = {r['sample'] for r in csv.DictReader(open(f'{BASE}/ebv_cn_1kGP_2504.tsv'), delimiter='\t') if r['sample'] in cn}
coinf = {s for s, v in wt.items() if v > 0.01}
clean_2502 = unrel_2502 - coinf

ir1, cp_d, uniq_d = {}, {}, {}
for line in open(f'{BASE}/metadata/all_ir1_2502.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 5 and f[4] != 'NA':
        try:
            ir1[f[0]] = float(f[4]) * RU
            cp_d[f[0]] = float(f[2])
            uniq_d[f[0]] = float(f[3])
        except ValueError: pass

# === 8. BCR CLONALITY: GERMLINE/PUBLIC PoN vs HARD CAP <= 20 ===
print('=== 8. WGS BCR CLONALITY: COHORT CDR3 PoN + V+J FILTER vs HARD READ-COUNT CAP ===')
# First consolidate or scan WGS BCR reports:
# Notice how fast it is if we check why >30 read rows failed RNA validation!
run2donor = {}
for m in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sequencing_library.metadata.MAGE.v1.0.txt'), delimiter='\t'):
    run2donor[m['SRA_accession']] = m['sample_kgpID']

rna_clones = defaultdict(dict)
for d in glob.glob(f'{BASE}/trust4/*/'):
    run = os.path.basename(d.rstrip('/'))
    rep = f'{d}{run}_report.tsv'
    if os.path.exists(rep) and run in run2donor:
        donor = run2donor[run]
        for i, line in enumerate(open(rep)):
            if i == 0: continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8: continue
            cnt, cdr3nt, cdr3aa, V, J = f[0], f[2], f[3], f[4], f[6]
            if J in ('.', '') or 'out_of_frame' in cdr3aa or '*' in cdr3aa or not cdr3aa or not cdr3nt or cdr3nt == '.': continue
            try: c = float(cnt)
            except ValueError: continue
            rna_clones[donor][cdr3nt] = rna_clones[donor].get(cdr3nt, 0.0) + c

# Load WGS clones for all donors that have RNA (plus a random 500 donor PoN for speed!)
wgs_files = glob.glob(f'{BASE}/bcr/*/*_report.tsv')
cdr3_donor_count = Counter()
wgs_raw = {}
for rep in wgs_files:
    s = os.path.basename(rep).replace('_report.tsv', '')
    rows = []
    seen_nt = set()
    for i, line in enumerate(open(rep)):
        if i == 0: continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8: continue
        cnt, cdr3nt, cdr3aa, V, D, J, C = f[0], f[2], f[3], f[4], f[5], f[6], f[7]
        if J in ('.', '') or 'out_of_frame' in cdr3aa or '*' in cdr3aa or not cdr3aa or not cdr3nt or cdr3nt == '.': continue
        try: c = float(cnt)
        except ValueError: continue
        has_V = V not in ('.', '')
        rows.append((c, cdr3nt, cdr3aa, has_V, V, J))
        seen_nt.add(cdr3nt)
    for nt in seen_nt:
        cdr3_donor_count[nt] += 1
    if rows:
        wgs_raw[s] = rows

shared_donors = sorted(set(rna_clones) & set(wgs_raw))
print(f'  Scanned {len(wgs_raw)} WGS BCR reports and {len(rna_clones)} RNA donors ({len(shared_donors)} paired)')

# Compare RNA recovery by read-count bin BEFORE vs AFTER requiring (has_V and cdr3_donor_count <= 2 [private to donor]):
for mode, filt_fn in [
    ('Old filter (J-only, no V check, no PoN)', lambda r: True),
    ('V+J required, no PoN', lambda r: r[3]),
    ('V+J + Private CDR3 PoN (<=2 donors in 2,502)', lambda r: r[3] and cdr3_donor_count[r[1]] <= 2)
]:
    bins = {'1-5': [0,0], '6-15': [0,0], '16-30': [0,0], '>30': [0,0]}
    for s in shared_donors:
        rset = rna_clones[s]
        # aggregate by cdr3nt
        by_nt = {}
        for c, nt, aa, has_V, V, J in wgs_raw[s]:
            if not filt_fn((c, nt, aa, has_V, V, J)): continue
            by_nt[nt] = by_nt.get(nt, 0.0) + c
        for nt, c in by_nt.items():
            b = '1-5' if c <= 5 else ('6-15' if c <= 15 else ('16-30' if c <= 30 else '>30'))
            bins[b][1] += 1
            if nt in rset: bins[b][0] += 1
    print(f'  {mode}:')
    for b in ('1-5', '6-15', '16-30', '>30'):
        k, n = bins[b]
        print(f'    count {b:<5}: {k:>4}/{n:<4} ({100*k/max(1,n):5.1f}% in RNA)')
