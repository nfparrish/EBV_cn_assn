"""IBS vs GBR EBV variant stratification -- the provenance test.

For each called site, compute the alternate-allele fraction per line from AD, then
compare IBS and GBR. Provenance predicts sites where one population carries a derived
allele the other lacks, co-segregating as a block.
"""
import gzip, statistics as st

V = '/work/nfp8/EBV/variants/pilot.vcf.gz'
pops = {}
for p in ('IBS', 'GBR'):
    for s in open(f'/work/nfp8/EBV/variants/{p}.ids').read().split():
        pops[s] = p

hdr, recs = None, []
with gzip.open(V, 'rt') as fh:
    for line in fh:
        if line.startswith('##'): continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#CHROM'):
            hdr = [c.split('/')[-1].replace('.chrEBV.bam', '') for c in f[9:]]
            continue
        recs.append(f)

fmt_cache = {}
print(f'{len(recs)} sites x {len(hdr)} samples')
print(f'IBS: {sum(1 for s in hdr if pops.get(s)=="IBS")}  GBR: {sum(1 for s in hdr if pops.get(s)=="GBR")}\n')
print(f'{"pos":>8} {"ref>alt":>10} {"IBS_AF":>7} {"GBR_AF":>7} {"diff":>6} {"IBS_n":>6} {"GBR_n":>6}')
rows = []
for f in recs:
    pos, ref, alt, fmt = f[1], f[3], f[4], f[8].split(':')
    try: adi = fmt.index('AD')
    except ValueError: continue
    af = {'IBS': [], 'GBR': []}
    for s, g in zip(hdr, f[9:]):
        p = pops.get(s)
        if not p: continue
        parts = g.split(':')
        if len(parts) <= adi or parts[adi] in ('.', ''): continue
        ad = [int(x) if x.isdigit() else 0 for x in parts[adi].split(',')]
        tot = sum(ad)
        if tot < 20: continue
        af[p].append(sum(ad[1:]) / tot)
    if len(af['IBS']) < 10 or len(af['GBR']) < 10: continue
    mi, mg = st.mean(af['IBS']), st.mean(af['GBR'])
    rows.append((abs(mi-mg), pos, f'{ref}>{alt}', mi, mg, len(af['IBS']), len(af['GBR'])))

rows.sort(reverse=True)
for d, pos, ra, mi, mg, ni, ng in rows:
    flag = '  <== STRATIFIED' if d > 0.3 else ''
    print(f'{pos:>8} {ra:>10} {mi:>7.3f} {mg:>7.3f} {d:>6.3f} {ni:>6} {ng:>6}{flag}')
print(f'\nsites with |IBS-GBR| allele-fraction difference > 0.3: {sum(1 for r in rows if r[0]>0.3)}')
print(f'sites analysed with adequate depth in both: {len(rows)}')
