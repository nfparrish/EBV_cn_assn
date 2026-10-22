"""Is the CN estimate inflated by tandem repeats collapsed in the reference?

Classical qPCR puts most LCLs under 10 episomes/cell (Raji 50-60, Namalwa 1-2), while our
depth-ratio estimate has median 33 and max 192. One candidate explanation: the CN numerator
excluded ONLY the B95-8 deletion (139,700-151,400). The repeat mask -- IR1, IR2, IR3, the
terminal repeats -- was derived later and applied only to VARIANT CALLING, so those regions
are still inside the CN numerator. Reads from the tandem copies of a repeat collapse onto the
single reference copy, inflating mean depth.
"""
import statistics as st
from collections import defaultdict

REG = [('IR1 (12.3-35.0 kb)', 12301, 35037),
       ('IR2', 38382, 39733),
       ('IR3 (EBNA1 Gly-Ala)', 95984, 96575),
       ('B95-8 deletion', 139724, 151551),
       ('terminal repeats', 169740, 171670)]

def region(p):
    for name, a, b in REG:
        if a <= p <= b: return name
    return 'unique sequence'

vals = defaultdict(list)
rows = []
with open('/work/nfp8/EBV/metadata/ebv_norm_depth.tsv') as fh:
    for line in fh:
        p, d = line.split('\t')
        p, d = int(p), float(d)
        rows.append((p, d)); vals[region(p)].append(d)

print('normalised depth by region (1.0 = per-sample mean):')
for name in [r[0] for r in REG] + ['unique sequence']:
    v = vals.get(name, [])
    if v: print(f'  {name:<24} mean {st.mean(v):6.3f}  median {st.median(v):6.3f}  over {len(v):6d} bp')

MASK_LO, MASK_HI = 139700, 151400          # what the CN numerator actually excluded
rep = [d for p, d in rows if not (MASK_LO <= p <= MASK_HI) and region(p) not in ('unique sequence', 'B95-8 deletion')]
uni = [d for p, d in rows if region(p) == 'unique sequence']
num = [d for p, d in rows if not (MASK_LO <= p <= MASK_HI)]
print(f'\nCN numerator as computed (excludes only the deletion):')
print(f'  repeat regions still inside: mean {st.mean(rep):.3f} over {len(rep)} bp')
print(f'  unique sequence            : mean {st.mean(uni):.3f} over {len(uni)} bp')
print(f'  numerator as computed      : mean {st.mean(num):.3f}')
print(f'\n  >>> inflation from including repeats: {st.mean(num)/st.mean(uni):.3f}x')
print(f'  >>> median CN 33 would become {33*st.mean(uni)/st.mean(num):.1f} on unique sequence only')
