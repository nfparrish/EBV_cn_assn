"""RNA vs WGS clonotype concordance, with the count-cap artifact filter.

KEY METHODOLOGICAL POINT: on genomic DNA, TRUST4 read count is INVERTED as a confidence signal.
Germline sequence sits at 2 copies x 30x and so attracts consistently HIGH counts, while a genuine
clonal rearrangement at 20% of cells yields only ~6 reads. Measured against RNA expression:
  1-5 reads 73.7% recovered | 6-15 84.6% | 16-30 54.8% | >30 only 1.0% (1/102).
So the filter must be an UPPER bound on count, not a lower bound -- the opposite of the usual
convention, and the opposite of what ranking by count implies.
"""
import csv, glob, os, random, statistics as st
from collections import defaultdict
CAP = 20

def load(path, cap=None):
    out = {}
    try:
        for i, line in enumerate(open(path)):
            if i == 0: continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8: continue
            if f[6] in ('.', '') or 'out_of_frame' in f[3] or '*' in f[3] or not f[3]: continue
            if not f[2] or f[2] == '.': continue
            try: c = float(f[0])
            except ValueError: continue
            if cap is not None and c > cap: continue
            out[f[2]] = out.get(f[2], 0.0) + c
    except FileNotFoundError: pass
    return out

r2d = {m['SRA_accession']: m['sample_kgpID'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/'
         'sequencing_library.metadata.MAGE.v1.0.txt'), delimiter='\t')}
rna = defaultdict(dict)
for d in glob.glob('/work/nfp8/EBV/trust4/*/'):
    run = os.path.basename(d.rstrip('/'))
    if run in r2d: rna[r2d[run]].update(load(f'{d}{run}_report.tsv'))
wgs = {}
for d in glob.glob('/work/nfp8/EBV/bcr/*/'):
    s = os.path.basename(d.rstrip('/'))
    c = load(f'{d}{s}_report.tsv', cap=CAP)
    if c: wgs[s] = c
shared = sorted(set(rna) & set(wgs))
print(f'paired donors: {len(shared)}   (WGS clones capped at <={CAP} reads)')
print(f'median WGS clones/donor after cap: {st.median([len(wgs[s]) for s in shared]):.0f}\n')

obs = [sum(1 for k in wgs[s] if k in rna[s])/len(wgs[s]) for s in shared]
random.seed(0)
ctl = []
for s in shared:
    o = random.choice([x for x in shared if x != s])
    ctl.append(sum(1 for k in wgs[s] if k in rna[o])/len(wgs[s]))
print('fraction of WGS clones recovered in RNA:')
print(f'  SAME donor  : mean {st.mean(obs):.3f}  median {st.median(obs):.3f}')
print(f'  CROSS donor : mean {st.mean(ctl):.3f}  median {st.median(ctl):.3f}')
print(f'  enrichment  : {st.mean(obs)/max(st.mean(ctl),1e-9):.0f}x')

t1 = t5 = inr = n = 0
for s in shared:
    w, r = wgs[s], rna[s]
    if not w or not r: continue
    n += 1
    wt = max(w, key=lambda k: w[k])
    order = sorted(r, key=lambda k: -r[k])
    if wt in r:
        inr += 1
        i = order.index(wt)
        t1 += (i == 0); t5 += (i < 5)
print(f'\ndominant-clone concordance after the cap (n={n}):')
print(f'  top WGS clone present in RNA      : {inr}/{n} ({100*inr/n:.0f}%)')
print(f'  top WGS clone is top RNA clone    : {t1}/{n} ({100*t1/n:.0f}%)')
print(f'  top WGS clone in RNA top 5        : {t5}/{n} ({100*t5/n:.0f}%)')

ranks = []
for s in shared:
    order = sorted(rna[s], key=lambda k: -rna[s][k])
    pos = {k: i/len(order) for i, k in enumerate(order)}
    ranks += [pos[k] for k in wgs[s] if k in pos]
print(f'\nRNA-abundance rank of recovered WGS clones: median {st.median(ranks):.3f} '
      f'(0 = most abundant, 0.5 = random; n={len(ranks)})')
