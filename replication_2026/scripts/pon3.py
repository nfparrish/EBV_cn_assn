"""Does derived-variant carriage structure by population across all 26 groups?

Provenance predicts population-structured viral variants; host-driven selection predicts
the same, so this cannot distinguish them on its own -- but its ABSENCE would refute both.
Permutation null: shuffle population labels, preserving group sizes and site frequencies.
"""
import csv, math, random, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'
wt, pop = {}, {}
with open(f'{BASE}/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try: wt[r['sample']] = float(r['wt_fraction'])
        except ValueError: pass
for m in csv.DictReader(open('/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/'
                             'sample_library_info/sample.metadata.MAGE.v1.0.txt'), delimiter='\t'):
    pop[m['sample_kgpID']] = m['population']
coinf = {s for s, v in wt.items() if v > 0.01}

site_car, samples = defaultdict(set), set()
with open(f'{BASE}/variants/all_calls.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) < 5: continue
        samples.add(p[0]); site_car[(int(p[1]), p[2], p[3])].add(p[0])
clean = samples - coinf

cohort = [s for s in clean if s in pop]
bypop = defaultdict(list)
for s in cohort: bypop[pop[s]].append(s)
bypop = {k: v for k, v in bypop.items() if len(v) >= 20}
cohort = [s for s in cohort if pop[s] in bypop]
NC = len(clean)
print(f'{len(cohort)} clean MAGE donors across {len(bypop)} populations\n')

derived = [(site, car & set(cohort)) for site, car in site_car.items()
           if 3 <= len(car & clean) <= 0.5*NC]
derived = [(s, c) for s, c in derived if len(c) >= 5]
print(f'derived sites with >=5 carriers in the MAGE subset: {len(derived)}')

def chi2(site_carriers):
    n = len(cohort); k = len(site_carriers)
    if k == 0 or k == n: return 0.0
    x = 0.0
    for p, members in bypop.items():
        obs = len(site_carriers & set(members)); exp = len(members) * k / n
        if exp > 0: x += (obs-exp)**2/exp
    return x

obs = [(chi2(c), site, len(c)) for site, c in derived]
obs.sort(reverse=True)
random.seed(0)
null = []
for _ in range(200):
    shuf = cohort[:]; random.shuffle(shuf)
    remap = {}; i = 0
    tmp = {}
    for p, members in bypop.items():
        tmp[p] = set(shuf[i:i+len(members)]); i += len(members)
    saved = dict(bypop); bypop.clear(); bypop.update({k: list(v) for k, v in tmp.items()})
    null.append(max(chi2(c) for _, c in derived))
    bypop.clear(); bypop.update(saved)

thr = sorted(null)[int(0.95*len(null))]
nsig = sum(1 for x, _, _ in obs if x > thr)
print(f'max-chi2 permutation threshold (95th pct of 200 label shuffles): {thr:.1f}')
print(f'sites exceeding it: {nsig}\n')
print(f'  {"site":<24} {"chi2":>8} {"carriers":>9}  top populations')
for x, site, nc in obs[:10]:
    car = dict(derived)[site]
    fr = sorted(((len(car & set(m))/len(m), p) for p, m in bypop.items()), reverse=True)[:3]
    flag = ' *' if x > thr else ''
    print(f'  chrEBV:{site[0]}{site[1]}>{site[2]:<9} {x:>8.1f} {nc:>9}  '
          + ', '.join(f'{p} {100*f:.0f}%' for f, p in fr) + flag)
