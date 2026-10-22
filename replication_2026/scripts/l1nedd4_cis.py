"""Cis test: does L1-NEDD4 genotype track NEDD4 expression in MAGE LCLs?

Kojima et al. (Nat Genet 2023) reported the cis effect in GTEx FIBROBLASTS, with exon 9
(short-isoform-specific) most strongly associated. LCLs are a different cell type, so this is
a replication attempt in a new context -- and a validation that genotype extraction, sample
matching and the expression tables all line up before anything is built on them.

Insertion: 1kgp_3202_v38_ins_2847, chr15:55,958,935, L1HS, SVLEN=6050, AC=1957/6404 (AF 0.306).
NEDD4 = ENSG00000069869, TSS chr15:55,993,659, minus strand.
"""
import gzip, math, statistics as st, csv
from collections import defaultdict

V = '/hpc/group/parrishlab/refs/megane_1kgp/1000GP.GRCh38_3202.ME_insertions.PASS.vcf.gz'
TARGET = '1kgp_3202_v38_ins_2847'
gt = {}
with gzip.open(V, 'rt') as fh:
    for line in fh:
        if line.startswith('##'): continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            samples = f[9:]; continue
        if f[2] != TARGET: continue
        for s, g in zip(samples, f[9:]):
            a = g.split(':')[0].replace('|', '/')
            if a in ('0/0',): gt[s] = 0
            elif a in ('0/1', '1/0'): gt[s] = 1
            elif a == '1/1': gt[s] = 2
        break
print(f'genotyped samples: {len(gt)}')
d = defaultdict(int)
for v in gt.values(): d[v] += 1
n = sum(d.values()); af = (d[1] + 2*d[2]) / (2*n)
print(f'  0/0={d[0]}  0/1={d[1]}  1/1={d[2]}   allele frequency={af:.3f}')

BED = ('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/QTL_results/eQTL_results/'
       'expression_quants/inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz')
with gzip.open(BED, 'rt') as fh:
    hdr = next(fh).rstrip('\n').split('\t'); ids = hdr[4:]
    row = None
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if f[3].startswith('ENSG00000069869'): row = f; break
expr = {s: float(v) for s, v in zip(ids, row[4:])}
shared = [s for s in ids if s in gt]
print(f'\nMAGE donors with both genotype and expression: {len(shared)}')

groups = defaultdict(list)
for s in shared: groups[gt[s]].append(expr[s])
print(f'\n{"genotype":<12} {"n":>5} {"mean inv-norm NEDD4":>22} {"sd":>7}')
for g in (0, 1, 2):
    v = groups[g]
    if v: print(f'{g:<12} {len(v):>5} {st.mean(v):>22.4f} {st.pstdev(v):>7.3f}')

# additive model: expression ~ dosage
x = [gt[s] for s in shared]; y = [expr[s] for s in shared]
mx, my = st.mean(x), st.mean(y)
dx = [a-mx for a in x]; dy = [b-my for b in y]
sxx = sum(a*a for a in dx)
beta = sum(a*b for a, b in zip(dx, dy)) / sxx
r = sum(a*b for a, b in zip(dx, dy)) / math.sqrt(sxx * sum(b*b for b in dy))
nn = len(shared)
t = abs(r) * math.sqrt((nn-2)/(1-r*r)) if abs(r) < 1 else 99
p = math.erfc(t/math.sqrt(2))
print(f'\nadditive dosage model (n={nn}):')
print(f'  beta = {beta:+.4f} inv-norm units per allele')
print(f'  r    = {r:+.4f}   p ~ {p:.3e}')
print(f'  direction: {"HIGHER NEDD4 in carriers" if beta>0 else "LOWER NEDD4 in carriers"}')
