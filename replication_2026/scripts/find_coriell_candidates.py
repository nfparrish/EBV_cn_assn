#!/usr/bin/env python3
from collections import defaultdict

ko = {}
with open('/work/nfp8/EBV/metadata/dennis_ko_lcls.tsv') as f:
    h = f.readline().strip().split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        while len(p)<len(h): p.append('')
        ko[p[0]] = dict(zip(h, p))

cn30, wt = {}, {}
with open('/work/nfp8/EBV/metadata/all_cn_3200.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=10: cn30[p[0]], wt[p[0]] = float(p[3]), float(p[9])

ir1 = {}
with open('/work/nfp8/EBV/metadata/all_ir1_2502.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6: ir1[p[0]] = float(p[5])*7.40

mage = {}
with open('/work/nfp8/EBV/metadata/all_mage_ebv_rna_779.tsv') as f:
    mh = f.readline().strip().split('\t')
    for line in f:
        p = line.strip().split('\t')
        mage[p[1]] = dict(zip(mh, p))

ped_pop = {}
with open('/work/nfp8/EBV/metadata/ped3202.txt') as f:
    for line in f:
        p = line.strip().split()
        if len(p)>=7 and p[0]!='FamilyID': ped_pop[p[1]] = p[5]

calls = defaultdict(list)
by_site = defaultdict(list)
with open('/work/nfp8/EBV/variants/all_calls.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6 and p[0]!='sample':
            s, pos, ref, alt, af, dp = p[0], int(p[1]), p[2], p[3], float(p[4]), int(p[5])
            calls[s].append((pos, ref, alt, af, dp))
            if wt.get(s, 0) <= 0.01:
                by_site[(pos, ref, alt)].append((af, dp, s))

print('=== A. CLEAN B95-8 HIGH-VAF (>=0.15) NON-REPEAT VARIANTS IN DENNIS KO 480 LCLs ===')
for k in sorted(ko):
    if wt.get(k, 0) > 0.01: continue
    for pos, ref, alt, af, dp in calls.get(k, []):
        if af >= 0.15 and not (7400 <= pos <= 8100) and not (12200 <= pos <= 35300) and not (36380 <= pos <= 36510) and not (38000 <= pos <= 41000) and not (95000 <= pos <= 98000) and not (166000 <= pos <= 172000):
            d = ko[k]
            m_ie = f"{float(mage[k]['lytic_ie_rpkm']):.1f}" if k in mage else '-'
            print(f"  {k} ({d['pop']}, set={d['dataset']}, LN2={d['ln2'] or '-'}) | {pos}:{ref}>{alt} | VAF={af:.3f} (dp={dp}) | CN={cn30.get(k,0):.1f} | IR1={ir1.get(k,0):.2f} | MAGE_IE={m_ie}")

print('\n=== B. TOP HIGHEST-VAF CLEAN B95-8 CARRIERS OF 12,129:T>A ACROSS ALL 2,449 CORIELL 1kGP ===')
for af, dp, s in sorted(by_site[(12129, 'T', 'A')], reverse=True)[:18]:
    in_ko = 'YES (' + ko[s]['dataset'] + ', LN2=' + (ko[s]['ln2'] or '-') + ')' if s in ko else 'Coriell order'
    m_ie = f"{float(mage[s]['lytic_ie_rpkm']):.1f}" if s in mage else '-'
    print(f"  {s} ({ped_pop.get(s,'?')}) | VAF={af:.3f} (dp={dp}) | CN={cn30.get(s,0):.1f} | IR1={ir1.get(s,0):.2f} | MAGE_IE={m_ie} | Ko={in_ko}")

print('\n=== C. ALL CLEAN B95-8 CARRIERS OF RARE Zp (90,944-91,500, excl 91,192) & Rp (92,896-93,500) PROMOTER VARIANTS ===')
for (pos, ref, alt), lst in sorted(by_site.items()):
    if len(lst) <= 10 and ((90944 <= pos <= 91500 and pos != 91192) or (92896 <= pos <= 93500)):
        reg = 'Zp_BZLF1p' if pos <= 91500 else 'Rp_BRLF1p'
        for af, dp, s in sorted(lst, reverse=True):
            in_ko = 'YES (' + ko[s]['dataset'] + ', LN2=' + (ko[s]['ln2'] or '-') + ')' if s in ko else 'Coriell'
            m_ie = f"{float(mage[s]['lytic_ie_rpkm']):.1f}" if s in mage else '-'
            print(f"  {reg} {pos}:{ref}>{alt} | {s} ({ped_pop.get(s,'?')}) | VAF={af:.3f} (dp={dp}) | CN={cn30.get(s,0):.1f} | IR1={ir1.get(s,0):.2f} | MAGE_IE={m_ie} | Ko={in_ko}")
