"""If a shared MIXED INOCULUM explains the provenance signal, carrier allele fractions
should be more similar WITHIN a population than between -- every line made from one stock
inherits that stock's variant frequency. Recurrent mutation predicts no such clustering.
"""
import csv, math, statistics as st
from collections import defaultdict

pop = {m['sample_kgpID']: m['population'] for m in csv.DictReader(
    open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/'
         'sample.metadata.MAGE.v1.0.txt'), delimiter='\t')}

af = defaultdict(dict)
with open('/work/nfp8/EBV/variants/all_calls.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 5 and p[1] in ('91192', '132965'):
            try: af[p[1]][p[0]] = float(p[4])
            except ValueError: pass

for site in ('91192', '132965'):
    d = af[site]
    g = defaultdict(list)
    for s, v in d.items():
        if s in pop: g[pop[s]].append(v)
    g = {k: v for k, v in g.items() if len(v) >= 3}
    allv = [v for vs in g.values() for v in vs]
    if len(allv) < 6: continue
    within = st.mean([st.pstdev(v) for v in g.values() if len(v) > 1])
    total = st.pstdev(allv)
    print(f'chrEBV:{site}  carriers in pops with n>=3: {len(allv)} across {len(g)} populations')
    print(f'  mean WITHIN-population sd of AF : {within:.4f}')
    print(f'  TOTAL sd of AF                  : {total:.4f}')
    print(f'  ratio within/total              : {within/total:.3f}   '
          f'({"clustered by population -> shared stock" if within/total < 0.8 else "not clustered"})')
    for k, v in sorted(g.items(), key=lambda x: -len(x[1]))[:5]:
        print(f'    {k:<5} n={len(v):<3} median AF {st.median(v):.3f}  sd {st.pstdev(v):.3f}')
    print()
