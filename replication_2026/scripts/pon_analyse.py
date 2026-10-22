"""Panel-of-normals analysis of EBV intrahost variants across 2,502 LCLs.

Every line descends from one B95-8 ancestor, so the cohort is its own panel of normals:
a site variable at similar allele fraction in a large fraction of samples is systematic
error (mapping, or B95-8-vs-NC_007605 reference difference); a site variable in a few
lines is a candidate DERIVED variant arising since transformation.

Co-infected lines (wt_fraction > 1%) carry a second divergent genome and would dominate
any derived-variant analysis, so they are partitioned out and reported separately.
"""
import csv, math, statistics as st
from collections import defaultdict

BASE = '/work/nfp8/EBV'

wt, cn, pop = {}, {}, {}
with open(f'{BASE}/ebv_cn_1kGP_2504.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            wt[r['sample']] = float(r['wt_fraction']); cn[r['sample']] = float(r['ebv_cn_depth_ratio'])
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
        s, pos, ref, alt, af = p[0], p[1], p[2], p[3], p[4]
        samples.add(s)
        try: site_af[(int(pos), ref, alt)][s] = float(af)
        except ValueError: pass

N = len(samples)
clean = samples - coinf
print(f'samples with calls: {N}   co-infected (wt>1%): {len(samples & coinf)}   clean: {len(clean)}')
print(f'distinct variant sites: {len(site_af)}\n')

rows = []
for site, d in site_af.items():
    nc = sum(1 for s in d if s in clean)
    afs = [v for s, v in d.items() if s in clean]
    rows.append((site, len(d), nc, st.median(afs) if afs else 0.0,
                 st.pstdev(afs) if len(afs) > 1 else 0.0,
                 sum(1 for s in d if s in coinf)))

NC = len(clean)
def band(nc):
    f = nc / NC
    if f > 0.50: return 'systematic (>50% of lines)'
    if f > 0.05: return 'common (5-50%)'
    if f > 0.01: return 'low-frequency (1-5%)'
    if nc >= 3:  return 'rare (3+ lines, <1%)'
    return 'private (1-2 lines)'

groups = defaultdict(list)
for r in rows: groups[band(r[2])].append(r)
order = ['systematic (>50% of lines)', 'common (5-50%)', 'low-frequency (1-5%)',
         'rare (3+ lines, <1%)', 'private (1-2 lines)']
print(f'{"class":<28} {"sites":>7} {"median AF":>10} {"median sd(AF)":>14}')
for k in order:
    g = groups.get(k, [])
    if not g: continue
    print(f'{k:<28} {len(g):>7} {st.median([x[3] for x in g]):>10.3f} '
          f'{st.median([x[4] for x in g]):>14.3f}')

print('\n--- the systematic sites (candidate reference/mapping artifacts) ---')
for site, n, nc, med, sd, nci in sorted(groups.get('systematic (>50% of lines)', []),
                                        key=lambda x: -x[2])[:15]:
    print(f'  chrEBV:{site[0]:<7} {site[1]}>{site[2]:<12} in {nc}/{NC} clean lines '
          f'({100*nc/NC:.0f}%)  median AF {med:.3f}  sd {sd:.3f}')

print('\n--- enrichment of calls in co-infected lines (sanity check) ---')
per_clean = sum(1 for s, d in ((s, d) for site, d in site_af.items() for s in d) if s in clean)
cc = sum(len([s for s in d if s in coinf]) for d in site_af.values())
nn = sum(len([s for s in d if s in clean]) for d in site_af.values())
if len(samples & coinf) and len(clean):
    print(f'  calls per co-infected line: {cc/len(samples & coinf):.1f}')
    print(f'  calls per clean line      : {nn/len(clean):.1f}')

# derived variants = polymorphic in clean lines, excluding systematic sites
derived = [r for r in rows if r[2] >= 3 and r[2] / NC <= 0.50]
print(f'\nderived candidate sites (3+ clean lines, <=50%): {len(derived)}')
with open(f'{BASE}/variants/pon_sites.tsv', 'w') as out:
    out.write('pos\tref\talt\tn_all\tn_clean\tfreq_clean\tmedian_af\tsd_af\tn_coinfected\tclass\n')
    for site, n, nc, med, sd, nci in sorted(rows, key=lambda x: -x[2]):
        out.write(f'{site[0]}\t{site[1]}\t{site[2]}\t{n}\t{nc}\t{nc/NC:.4f}\t'
                  f'{med:.4f}\t{sd:.4f}\t{nci}\t{band(nc)}\n')
print(f'wrote {BASE}/variants/pon_sites.tsv')
