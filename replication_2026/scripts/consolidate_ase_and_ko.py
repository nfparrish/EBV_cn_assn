#!/usr/bin/env python3
import glob
from collections import defaultdict

run_to_donor = {}
with open('/work/nfp8/EBV/metadata/all_mage_ebv_rna_779.tsv') as f:
    h = f.readline().strip().split('\t')
    for l in f:
        p = l.strip().split('\t')
        run_to_donor[p[0]] = p[1]

wgs_calls = defaultdict(dict)
with open('/work/nfp8/EBV/variants/all_calls.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=6 and p[0]!='sample':
            wgs_calls[p[0]][(int(p[1]), p[2], p[3])] = (float(p[4]), int(p[5]))

with open('/work/nfp8/EBV/metadata/all_mage_viral_ase_59.tsv', 'w') as out:
    out.write('run\tdonor\tpos\tref\talt\tdna_vaf\tdna_dp\trna_ref\trna_alt\trna_vaf\n')
    for fn in sorted(glob.glob('/work/nfp8/EBV/ase_mage/*.ase.tsv')):
        with open(fn) as f:
            f.readline()
            for l in f:
                run, pos, ref, alt, r_ref, r_alt = l.strip().split('\t')
                pos, r_ref, r_alt = int(pos), int(r_ref), int(r_alt)
                donor = run_to_donor.get(run, '')
                if (pos, ref, alt) in wgs_calls[donor] or r_alt > 0:
                    daf, ddp = wgs_calls[donor].get((pos, ref, alt), (0.0, 0))
                    rtot = r_ref + r_alt
                    rvaf = r_alt / rtot if rtot > 0 else 0.0
                    if daf > 0 or (rtot >= 10 and r_alt >= 2):
                        out.write(f'{run}\t{donor}\t{pos}\t{ref}\t{alt}\t{daf:.4f}\t{ddp}\t{r_ref}\t{r_alt}\t{rvaf:.4f}\n')

with open('/work/nfp8/EBV/metadata/all_mage_viral_splices_59.tsv', 'w') as out:
    out.write('run\tdonor\tsplice_class\tcount\n')
    for fn in sorted(glob.glob('/work/nfp8/EBV/ase_mage/*.spl.tsv')):
        with open(fn) as f:
            f.readline()
            for l in f:
                run, sc, cnt = l.strip().split('\t')
                d = run_to_donor.get(run, '')
                out.write(f'{run}\t{d}\t{sc}\t{cnt}\n')

# Also generate annotated 480-row Dennis Ko LCL master table
cn30, wt = {}, {}
with open('/work/nfp8/EBV/metadata/all_cn_3200.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=10: cn30[p[0]], wt[p[0]] = float(p[3]), float(p[9])

ont_cn = {}
with open('/work/nfp8/EBV/metadata/all_ont_904.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=6: ont_cn[p[0]] = float(p[5])

low_cn = {}
with open('/work/nfp8/EBV/metadata/all_lowcov_2499.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=6: low_cn[p[0]] = float(p[5])

ir1_30 = {}
with open('/work/nfp8/EBV/metadata/all_ir1_2502.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=6: ir1_30[p[0]] = float(p[5])*7.40

ir1_ont = {}
with open('/work/nfp8/EBV/metadata/all_ont_ir1_907.tsv') as f:
    for l in f:
        p = l.strip().split('\t')
        if len(p)>=4: ir1_ont[p[0]] = float(p[3])

mage = {}
with open('/work/nfp8/EBV/metadata/all_mage_ebv_rna_779.tsv') as f:
    mh = f.readline().strip().split('\t')
    for l in f:
        p = l.strip().split('\t')
        mage[p[1]] = dict(zip(mh, p))

with open('/work/nfp8/EBV/metadata/dennis_ko_lcls.tsv') as f, open('/work/nfp8/EBV/metadata/dennis_ko_lcls_annotated_480.tsv', 'w') as out:
    h = f.readline().strip().split('\t')
    out.write('\t'.join(h + ['cn_30x', 'cn_ont', 'cn_lowcov', 'wt_fraction', 'ir1_30x', 'ir1_ont', 'mage_lytic_ie_rpkm', 'mage_latent_rpkm', 'mage_bzlf1_spliced_reads', 'key_viral_variants']) + '\n')
    for l in f:
        p = l.rstrip('\n').split('\t')
        while len(p)<len(h): p.append('')
        k = p[0]
        c30 = f'{cn30[k]:.2f}' if k in cn30 else ''
        cont = f'{ont_cn[k]:.2f}' if k in ont_cn else ''
        clow = f'{low_cn[k]:.2f}' if k in low_cn else ''
        wf = f'{wt[k]:.4f}' if k in wt else ''
        i30 = f'{ir1_30[k]:.2f}' if k in ir1_30 else ''
        iont = f'{ir1_ont[k]:.2f}' if k in ir1_ont else ''
        mie = f"{float(mage[k]['lytic_ie_rpkm']):.2f}" if k in mage else ''
        mlat = f"{float(mage[k]['latent_mrna_rpkm']):.2f}" if k in mage else ''
        mz = mage[k]['bzlf1_spliced_reads'] if k in mage else ''
        kv = []
        for (pos, ref, alt), (af, dp) in sorted(wgs_calls[k].items()):
            if pos in (12129, 163858, 91192, 132965, 168963, 36470, 36494) or (90944<=pos<=91500) or (92896<=pos<=93500):
                kv.append(f'{pos}:{ref}>{alt}({af:.3f})')
        out.write('\t'.join(p + [c30, cont, clow, wf, i30, iont, mie, mlat, mz, ';'.join(kv)]) + '\n')
print('Consolidated successfully.')
