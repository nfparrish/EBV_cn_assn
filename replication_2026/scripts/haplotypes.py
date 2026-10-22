"""Reconstruct EBV stock haplotypes from variant co-segregation across 2,449 clean LCLs.

All lines descend from B95-8, so derived variants that co-segregate across many
independently established lines cannot be independent recurrent mutations -- they were
inherited together from a shared inoculum. Clustering sites by co-occurrence therefore
reconstructs the variant content of the virus stocks used to make the panels.

Co-infected lines (wt>1%) are excluded: they share wild-type variants and would form one
large spurious cluster.
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

car, samples = defaultdict(set), set()
with open(f'{BASE}/variants/all_calls.tsv') as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) < 5: continue
        samples.add(p[0]); car[(int(p[1]), p[2], p[3])].add(p[0])
clean = samples - coinf
N = len(clean)

MINCAR, MAXFRAC = 5, 0.50
sites = {s: c & clean for s, c in car.items()}
sites = {s: c for s, c in sites.items() if MINCAR <= len(c) <= MAXFRAC * N}
keys = sorted(sites)
print(f'{N} clean lines; {len(keys)} derived sites with >={MINCAR} carriers\n')

def logp_hyper(k, n1, n2, n):
    """-log10 P(overlap >= k) under independence, via a normal approximation."""
    mu = n1 * n2 / n
    if k <= mu: return 0.0
    var = mu * (1 - n1/n) * (1 - n2/n) * n/(n-1) if n > 1 else 1
    if var <= 0: return 0.0
    z = (k - mu - 0.5) / math.sqrt(var)
    return -math.log10(max(math.erfc(z / math.sqrt(2)) / 2, 1e-300))

# pairwise co-segregation, Bonferroni over all tested pairs
npairs = len(keys) * (len(keys) - 1) / 2
THR = -math.log10(0.01 / max(npairs, 1))
edges = defaultdict(set)
pairinfo = []
for i, a in enumerate(keys):
    ca = sites[a]
    for b in keys[i+1:]:
        cb = sites[b]
        k = len(ca & cb)
        if k < 4: continue
        lp = logp_hyper(k, len(ca), len(cb), N)
        if lp > THR:
            j = k / len(ca | cb)
            edges[a].add(b); edges[b].add(a)
            pairinfo.append((lp, j, a, b, k, len(ca), len(cb)))
print(f'Bonferroni threshold -log10p > {THR:.1f}; significant co-segregating pairs: {len(pairinfo)}')

# connected components = candidate stock haplotypes
seen, comps = set(), []
for s in keys:
    if s in seen or s not in edges: continue
    stack, comp = [s], set()
    while stack:
        x = stack.pop()
        if x in comp: continue
        comp.add(x); seen.add(x)
        stack.extend(edges[x] - comp)
    if len(comp) >= 2: comps.append(comp)
comps.sort(key=len, reverse=True)
print(f'haplotype clusters (>=2 sites): {len(comps)}\n')

mage = {s for s in clean if s in pop}
bypop = defaultdict(set)
for s in mage: bypop[pop[s]].add(s)
bypop = {k: v for k, v in bypop.items() if len(v) >= 20}

print(f'{"#":<3} {"sites":>5} {"core lines":>10}  populations (carrier % of that population)')
for i, comp in enumerate(comps[:14], 1):
    # lines carrying a majority of the cluster's sites = cluster carriers
    cnt = defaultdict(int)
    for s in comp:
        for x in sites[s]: cnt[x] += 1
    core = {x for x, c in cnt.items() if c >= max(2, len(comp)/2)}
    fr = []
    for p, mem in bypop.items():
        f = len(core & mem) / len(mem)
        if f > 0.05: fr.append((f, p, len(core & mem), len(mem)))
    fr.sort(reverse=True)
    desc = ', '.join(f'{p} {100*f:.0f}% ({k}/{n})' for f, p, k, n in fr[:4]) or 'no population >5%'
    print(f'{i:<3} {len(comp):>5} {len(core):>10}  {desc}')
    for s in sorted(comp)[:6]:
        print(f'       chrEBV:{s[0]} {s[1]}>{s[2]} ({len(sites[s])} carriers)')
    if len(comp) > 6: print(f'       ... +{len(comp)-6} more sites')
    print()
