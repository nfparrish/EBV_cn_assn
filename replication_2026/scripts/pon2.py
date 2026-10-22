"""Panel-of-normals classification + population-structure test of DERIVED EBV variants.

Cohort-as-PoN: all 2,502 lines descend from one B95-8 ancestor, so a site variable at
similar AF in most lines is systematic error; one variable in a few lines is a candidate
variant derived since transformation.

Co-infected lines (wt_fraction>1%) carry a second divergent genome and are partitioned out.
NOTE: sites seen ONLY in co-infected lines must be their own class -- lumping them with
"private" (as a first pass did) misreports them, since they have zero clean carriers.
"""
import csv, math, statistics as st
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
site_af, samples = defaultdict(dict), set()
with open(f'{BASE}/variants/all_calls.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) < 5: continue
        samples.add(p[0])
        try: site_af[(int(p[1]), p[2], p[3])][p[0]] = float(p[4])
        except ValueError: pass
clean = samples - coinf
NC = len(clean)

def band(nc):
    if nc == 0:  return 'co-infected lines only'
    f = nc / NC
    if f > 0.50: return 'systematic (>50%)'
    if f > 0.05: return 'common (5-50%)'
    if f > 0.01: return 'low-frequency (1-5%)'
    if nc >= 3:  return 'rare (3+ lines)'
    return 'private (1-2 lines)'

groups = defaultdict(list)
for site, d in site_af.items():
    afs = [v for s, v in d.items() if s in clean]
    groups[band(len(afs))].append((site, len(afs), afs))
print(f'clean lines {NC}, co-infected {len(samples & coinf)}, distinct sites {len(site_af)}\n')
print(f'{"class":<26} {"sites":>7} {"median AF":>10}')
for k in ['systematic (>50%)','common (5-50%)','low-frequency (1-5%)','rare (3+ lines)',
          'private (1-2 lines)','co-infected lines only']:
    g = groups.get(k, [])
    if g:
        allaf = [a for _, _, afs in g for a in afs]
        print(f'{k:<26} {len(g):>7} {st.median(allaf) if allaf else float("nan"):>10.3f}')

# ---- where do the systematic sites sit? (mask-adjacency check)
MASK = [(12301,35037),(38382,39733),(95984,96575),(139724,151551),(169740,170650),(170764,171670)]
def near_mask(p, pad=2000):
    return any(a-pad <= p <= b+pad for a, b in MASK)
sysf = groups.get('systematic (>50%)', [])
print(f'\nsystematic sites within 2 kb of the empirical repeat mask: '
      f'{sum(1 for s,_,_ in sysf if near_mask(s[0]))}/{len(sysf)}')
hot = defaultdict(int)
for s, _, _ in sysf: hot[s[0]//1000] += 1
print('  positional clustering (kb bin: n sites): ' +
      ', '.join(f'{k}kb:{v}' for k, v in sorted(hot.items(), key=lambda x: -x[1])[:6]))

# ---- population structure of derived variants, IBS vs GBR and cohort-wide
derived = [(s, afs) for k in ('common (5-50%)','low-frequency (1-5%)','rare (3+ lines)')
           for s, n, afs in groups.get(k, [])]
print(f'\nderived candidate sites (3+ clean lines, <=50%): {len(derived)}')

carriers = {site: {s for s in site_af[site] if s in clean} for site, _ in derived}
bypop = defaultdict(set)
for s in clean:
    if s in pop: bypop[pop[s]].add(s)
bypop = {k: v for k, v in bypop.items() if len(v) >= 20}
print(f'populations with >=20 clean MAGE donors: {len(bypop)}'
      f' (total {sum(len(v) for v in bypop.values())} donors)')

def fisher_like(a, b, c, d):      # normal approx to a 2x2 log-odds test
    if min(a, b, c, d) == 0: a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
    lor = math.log((a*d)/(b*c)); se = math.sqrt(1/a+1/b+1/c+1/d)
    return lor, math.erfc(abs(lor/se)/math.sqrt(2))

if 'IBS' in bypop and 'GBR' in bypop:
    I, G = bypop['IBS'], bypop['GBR']
    print(f'\n--- IBS (n={len(I)}) vs GBR (n={len(G)}): derived-variant carrier frequency ---')
    res = []
    for site, _ in derived:
        ci = len(carriers[site] & I); cg = len(carriers[site] & G)
        if ci + cg == 0: continue
        lor, p = fisher_like(ci, len(I)-ci, cg, len(G)-cg)
        res.append((p, site, ci, cg, lor))
    res.sort()
    print(f'  sites tested: {len(res)}   nominal p<0.05: {sum(1 for r in res if r[0]<0.05)}'
          f'   expected by chance: {0.05*len(res):.1f}')
    print(f'  {"site":<22} {"IBS":>8} {"GBR":>8} {"log-OR":>8} {"p":>9} {"BH q":>9}')
    m = len(res)
    for k, (p, site, ci, cg, lor) in enumerate(res[:10], 1):
        q = min(1, p*m/k)
        print(f'  chrEBV:{site[0]}{site[1]}>{site[2]:<8} {ci:>3}/{len(I):<4} {cg:>3}/{len(G):<4} '
              f'{lor:>+8.2f} {p:>9.2e} {q:>9.2e}')
    sig = [r for k, r in enumerate(res, 1) if min(1, r[0]*m/k) < 0.05]
    print(f'\n  sites passing BH q<0.05: {len(sig)}')
