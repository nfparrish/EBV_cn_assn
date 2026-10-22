#!/usr/bin/env python3
from collections import defaultdict
import math

ko = {}
with open('/work/nfp8/EBV/metadata/dennis_ko_lcls.tsv') as f:
    h = f.readline().strip().split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        while len(p) < len(h): p.append('')
        d = dict(zip(h, p))
        ko[d['iid']] = d

cn30 = {}
wt_frac = {}
with open('/work/nfp8/EBV/metadata/all_cn_3200.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=10:
            cn30[p[0]] = float(p[3])
            wt_frac[p[0]] = float(p[9])

ont_cn = {}
with open('/work/nfp8/EBV/metadata/all_ont_904.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6: ont_cn[p[0]] = float(p[5])

low_cn = {}
with open('/work/nfp8/EBV/metadata/all_lowcov_2499.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6: low_cn[p[0]] = float(p[5])

ir1_30 = {}
with open('/work/nfp8/EBV/metadata/all_ir1_2502.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6: ir1_30[p[0]] = float(p[5]) * 7.40

ir1_ont = {}
with open('/work/nfp8/EBV/metadata/all_ont_ir1_907.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=4: ir1_ont[p[0]] = float(p[3])

mage = {}
with open('/work/nfp8/EBV/metadata/all_mage_ebv_rna_779.tsv') as f:
    mh = f.readline().strip().split('\t')
    for line in f:
        p = line.strip().split('\t')
        d = dict(zip(mh, p))
        mage[d['donor']] = d

calls = defaultdict(list)
with open('/work/nfp8/EBV/variants/all_calls.tsv') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p)>=6 and p[0] != 'sample':
            calls[p[0]].append((int(p[1]), p[2], p[3], float(p[4]), int(p[5])))

print('Dennis Ko 480 LCLs coverage across our datasets:')
print('  Total Dennis Ko LCLs:', len(ko))
print('  In 30x CN (3,200):', sum(1 for k in ko if k in cn30))
print('  In 30x IR1 / LoFreq Variants (2,502):', sum(1 for k in ko if k in ir1_30))
print('  In ONT CN / IR1 (907):', sum(1 for k in ko if k in ont_cn))
print('  In MAGE RNA-seq (730):', sum(1 for k in ko if k in mage))

def fmt_donor(k):
    d = ko[k]
    ln2 = d['ln2'] if d['ln2'] else '-'
    mem = d['member'] if d['member'] else 'unrel'
    m_ie = f"{float(mage[k]['lytic_ie_rpkm']):.1f}" if k in mage else '-'
    m_lat = f"{float(mage[k]['latent_mrna_rpkm']):.1f}" if k in mage else '-'
    m_z = mage[k]['bzlf1_spliced_reads'] if k in mage else '-'
    c30 = f"{cn30[k]:.1f}" if k in cn30 else '-'
    cont = f"{ont_cn[k]:.1f}" if k in ont_cn else '-'
    i30 = f"{ir1_30[k]:.2f}" if k in ir1_30 else '-'
    iont = f"{ir1_ont[k]:.2f}" if k in ir1_ont else '-'
    return f"{k} | pop={d['pop']} | sex={d['sex']} | fam={d['fam']}({mem}) | LN2={ln2} | set={d['dataset']} | CN(30x/ONT)={c30}/{cont} | IR1(30x/ONT)={i30}/{iont} | MAGE(IE/Lat/Zspl)={m_ie}/{m_lat}/{m_z}"

print('\n=== 1. DENNIS KO LCLs CARRYING chrEBV:12,129 T>A (5p W0/Cp-Wp Leader CN-QTL & Lytic IE eQTL) ===')
for k in sorted(ko):
    for pos, ref, alt, af, dp in calls.get(k, []):
        if pos == 12129 and alt == 'A':
            print(f"  VAF_DNA={af:.3f} (dp={dp:4d}) | {fmt_donor(k)}")

print('\n=== 2. DENNIS KO LCLs CARRYING chrEBV:163,858 G>T (BALF2 ssDNA-Binding Protein Cis/Trans Lytic Suppressor QTL) ===')
for k in sorted(ko):
    for pos, ref, alt, af, dp in calls.get(k, []):
        if pos == 163858 and alt == 'T':
            print(f"  VAF_DNA={af:.3f} (dp={dp:4d}) | {fmt_donor(k)}")

print('\n=== 3. DENNIS KO LCLs CARRYING Zp (90,944-91,500), Rp (92,896-93,500), OR oriP (7,315-9,312) VARIANTS ===')
for k in sorted(ko):
    for pos, ref, alt, af, dp in calls.get(k, []):
        reg = None
        if 90944 <= pos <= 91500: reg = 'Zp/BRLF1'
        elif 92896 <= pos <= 93500: reg = 'Rp_BRLF1p'
        elif 7315 <= pos <= 9312: reg = 'oriP'
        if reg:
            print(f"  {reg:<9s} {pos}:{ref}>{alt:<3s} VAF={af:.3f} (dp={dp:4d}) | {fmt_donor(k)}")

print('\n=== 4. DENNIS KO LCLs CARRYING CIS-ASE VARIANTS (168,963:T>G in LMP-1; 36,470:T>C or 36,494:G>C in EBNA-2; 132,965:G>GC in BXRF1) ===')
for k in sorted(ko):
    for pos, ref, alt, af, dp in calls.get(k, []):
        if pos in (168963, 36470, 36494) or (pos == 132965 and af >= 0.20):
            print(f"  {pos}:{ref}>{alt:<3s} VAF={af:.3f} (dp={dp:4d}) | {fmt_donor(k)}")

print('\n=== 5. DENNIS KO EXTREME IR1 REPEAT COPY NUMBER PAIRS WITHIN THE SAME POPULATION (Structural QTL Test) ===')
by_pop = defaultdict(list)
for k in ko:
    if k in ir1_30 and k in cn30 and wt_frac.get(k,0) <= 0.01:
        by_pop[ko[k]['pop']].append(k)
for pop, ks in sorted(by_pop.items()):
    ks_s = sorted(ks, key=lambda x: ir1_30[x])
    print(f"  --- Population {pop} (n={len(ks_s)} in 2,502 WGS) ---")
    for k in ks_s[:3] + ks_s[-3:]:
        print(f"    {fmt_donor(k)}")

print('\n=== 6. DENNIS KO SAME-FAMILY TRIO / SPOUSE PAIRS WITH MAXIMUM EBV CN OR IR1 DISCORDANCE ===')
by_fam = defaultdict(list)
for k, d in ko.items():
    if d['fam'] and not d['fam'].startswith('HG') and k in cn30:
        by_fam[d['fam']].append(k)
fam_diffs = []
for fam, ks in by_fam.items():
    if len(ks) >= 2:
        cns = [cn30[k] for k in ks]
        ratio = max(cns) / max(1.0, min(cns))
        ir_d = 0.0
        irs = [ir1_30[k] for k in ks if k in ir1_30]
        if len(irs) >= 2: ir_d = max(irs) - min(irs)
        fam_diffs.append((ratio, ir_d, fam, ks))
for ratio, ir_d, fam, ks in sorted(fam_diffs, reverse=True)[:10]:
    print(f"  Family {fam} (max/min CN ratio = {ratio:.1f}x, ΔIR1 = {ir_d:.2f}u):")
    for k in ks:
        vstr = ','.join(f'{pos}:{ref}>{alt}({af:.2f})' for pos,ref,alt,af,dp in calls.get(k,[]) if pos in (12129,163858,91192,132965,168963,36470,36494) or 90944<=pos<=91500 or 92896<=pos<=93500)
        print(f"    {fmt_donor(k)} | key_vars={vstr or 'none'}")
