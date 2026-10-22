import csv, glob, math, os, statistics as st
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

def rank(vals):
    o = sorted(range(len(vals)), key=lambda i: vals[i]); r = [0.0]*len(vals); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and vals[o[j+1]] == vals[o[i]]: j += 1
        for k in range(i, j + 1): r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r

def pearson(a, b):
    ma, mb = st.mean(a), st.mean(b)
    da = [x - ma for x in a]; db = [y - mb for y in b]
    den = math.sqrt(sum(x*x for x in da) * sum(y*y for y in db))
    return sum(x*y for x, y in zip(da, db)) / den if den else 0.0

def spearman(a, b):
    return pearson(rank(a), rank(b))

def pv(r, n):
    if n < 5 or abs(r) >= 1: return 0.0
    t = abs(r) * math.sqrt((n - 2) / (1 - r*r))
    return math.erfc(t / math.sqrt(2))

# 1. PoN-cleaned WGS BCR clonality vs CN and IR1 (NO count<=20 cap!)
wgs_files = glob.glob(f'{BASE}/bcr/*/*_report.tsv')
cdr3_donor_count = Counter()
wgs_raw = {}
for rep in wgs_files:
    s = os.path.basename(rep).replace('_report.tsv', '')
    by_nt = {}
    for i, line in enumerate(open(rep)):
        if i == 0: continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 8: continue
        cnt, cdr3nt, cdr3aa, V, J = f[0], f[2], f[3], f[4], f[6]
        if V in ('.', '') or J in ('.', '') or 'out_of_frame' in cdr3aa or '*' in cdr3aa or not cdr3aa or not cdr3nt or cdr3nt == '.':
            continue
        try: c = float(cnt)
        except ValueError: continue
        if c > 0:
            by_nt[cdr3nt] = by_nt.get(cdr3nt, 0.0) + c
    for nt in by_nt:
        cdr3_donor_count[nt] += 1
    if by_nt:
        wgs_raw[s] = by_nt

bcr_pon = {}
for s, by_nt in wgs_raw.items():
    cl = [c for nt, c in by_nt.items() if cdr3_donor_count[nt] <= 2]
    if len(cl) >= 2:
        tot = sum(cl); top = max(cl) / tot
        H = -sum((c/tot)*math.log(c/tot) for c in cl)
        bcr_pon[s] = (len(cl), top, H / math.log(len(cl)), max(cl))

ss_bcr = sorted(s for s in bcr_pon if s in clean_2502)
print(f'=== 8b. PoN-CLEANED WGS BCR CLONALITY (NO READ-COUNT CAP) vs CN & IR1 (n={len(ss_bcr)} clean lines) ===')
for idx, name in [(0, 'n_clones'), (1, 'top_clone_frac'), (2, 'norm_Shannon'), (3, 'max_clone_reads')]:
    v = [bcr_pon[s][idx] for s in ss_bcr]
    y_cn = [math.log2(cn[s]) for s in ss_bcr]
    r_cn = spearman(v, y_cn)
    y_ir = [ir1[s] for s in ss_bcr if s in ir1]
    v_ir = [bcr_pon[s][idx] for s in ss_bcr if s in ir1]
    r_ir = spearman(v_ir, y_ir)
    print(f'  {name:<16}: median={st.median(v):.3f} | vs log2(CN): rho={r_cn:+.4f} (p={pv(r_cn, len(ss_bcr)):.2e}) | vs IR1: rho={r_ir:+.4f} (p={pv(r_ir, len(v_ir)):.2e})')

# Also check RNA vs WGS dominant clone concordance with PoN (no cap<=20!)
run2donor = {m['SRA_accession']: m['sample_kgpID'] for m in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sequencing_library.metadata.MAGE.v1.0.txt'), delimiter='\t')}
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
            cnt, cdr3nt, cdr3aa, J = f[0], f[2], f[3], f[6]
            if J in ('.', '') or 'out_of_frame' in cdr3aa or '*' in cdr3aa or not cdr3aa or not cdr3nt or cdr3nt == '.': continue
            try: c = float(cnt)
            except ValueError: continue
            rna_clones[donor][cdr3nt] = rna_clones[donor].get(cdr3nt, 0.0) + c

shared_d = sorted(s for s in rna_clones if s in wgs_raw)
top1 = top5 = inrna = n_eval = 0
obs_f = []
for s in shared_d:
    w = {nt: c for nt, c in wgs_raw[s].items() if cdr3_donor_count[nt] <= 2}
    r = rna_clones[s]
    if not w or not r: continue
    n_eval += 1
    obs_f.append(sum(1 for k in w if k in r) / len(w))
    wtop = max(w, key=lambda k: w[k])
    rorder = sorted(r, key=lambda k: -r[k])
    if wtop in r:
        inrna += 1
        idx = rorder.index(wtop)
        if idx == 0: top1 += 1
        if idx < 5: top5 += 1
print(f'  PoN-cleaned WGS vs RNA concordance (n={n_eval} donors, NO count<=20 cap):')
print(f'    Mean WGS clones in RNA: {100*st.mean(obs_f):.1f}% (median {100*st.median(obs_f):.1f}%)')
print(f'    Top WGS clone in RNA  : {inrna}/{n_eval} ({100*inrna/n_eval:.1f}%), in RNA top 5: {top5}/{n_eval} ({100*top5/n_eval:.1f}%), IS top RNA clone: {top1}/{n_eval} ({100*top1/n_eval:.1f}%)')
