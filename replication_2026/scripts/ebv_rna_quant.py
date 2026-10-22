#!/usr/bin/env python3
"""Virologically rigorous EBV RNA-seq quantification for MAGE & GEUVADIS LCLs.

Solves three EBV-specific transcriptomic traps:
  1. Competitive host filtering: discards any candidate viral read whose alignment match
     to a human GENCODE v45 transcript (e.g., IL10 vs viral BCRF1, CSF1R vs BARF1, BCL2 vs
     BHRF1, human Gly-Ala repeats vs EBNA1 IR3) >= 90% of its EBV alignment span.
  2. Non-overlapping CDS windows: prevents polycistronic/intronic spans of EBNA-1..3C
     (11.3-97.6 kb) and BARTs from swallowing or colliding with lytic ORFs.
  3. Separate Pol III (EBER1/2) vs Pol II Latent mRNA vs IE vs Early vs Late Structural
     quantification, plus direct CIGAR 'N' splice-junction counting (separating spliced
     BZLF1 Zta mRNA from unspliced BRLF1 bicistronic read-through, and counting Cp/Wp/Qp
     EBNA splices).
"""
import csv, math, os, re, subprocess, sys

run_id, sam_path, paf_path, out_tsv, out_genes = sys.argv[1:6]

# Load library metadata (numReads, donor, batch, pop, RIN)
meta_path = '/work/nfp8/MAGE/metadata/zenodo/MAGE.v1.0.data/sample_library_info/sequencing_library.metadata.MAGE.v1.0.txt'
donor, batch, pop, cont, num_reads, rin = 'NA', 'NA', 'NA', 'NA', 25000000.0, 'NA'
if os.path.exists(meta_path):
    for r in csv.DictReader(open(meta_path), delimiter='\t'):
        if r['SRA_accession'] == run_id:
            donor = r['sample_kgpID']
            batch = r['batch']
            pop   = r['population']
            cont  = r['continentalGroup']
            try: num_reads = float(r['numReads']) * 2.0 # numReads is read pairs; total reads = 2 * numReads
            except ValueError: pass
            rin   = r['RIN']
            break

# 1. Parse host competitive PAF: record best human transcript match bases per qname
host_best = {}
if os.path.exists(paf_path):
    for line in open(paf_path):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 11:
            qname = f[0]
            mlen  = int(f[9]) # number of residue matches on human transcript
            if mlen > host_best.get(qname, 0):
                host_best[qname] = mlen

# 2. Curated non-overlapping EBV kinetic intervals (1-based inclusive on NC_007605.1 / chrEBV)
# Format: (gene_name, kinetic_class, start, end, strand)
# Kinetic classes:
#   LATENT_MRNA : strict latent Pol II protein-coding exons (excluding IR1 repeat & IR3 Gly-Ala repeat)
#   LATENT_BART : BART lncRNA exons outside B95-8 deletion
#   POL3_EBER   : Pol III non-polyadenylated EBER1 & EBER2 (excluded from mRNA ratios!)
#   LYTIC_IE    : Immediate-early lytic (BZLF1, BRLF1)
#   LYTIC_EARLY : Early lytic replication/regulation
#   LYTIC_LATE  : True late structural (capsid, tegument, envelope glycoproteins)
#   B958_DEL    : Inside B95-8 deletion (139,724-151,400: LF3, LF2, LF1) -> wild-type EBV RNA detector
GENES = [
    # LATENT_MRNA
    ('LMP-2B_ex',   'LATENT_MRNA',     59,   1680, '+'),
    ('EBNA-2',      'LATENT_MRNA',  36216,  37679, '+'),
    ('EBNA-3A',     'LATENT_MRNA',  79955,  82877, '+'),
    ('EBNA-3B_3C',  'LATENT_MRNA',  83065,  89135, '+'),
    ('EBNA-1_N',    'LATENT_MRNA',  95662,  95950, '+'), # 5' of IR3 Gly-Ala repeat
    ('EBNA-1_C',    'LATENT_MRNA',  96600,  97587, '+'), # 3' of IR3 Gly-Ala repeat
    ('LMP-2A_ex1',  'LATENT_MRNA', 166103, 166458, '+'),
    ('LMP-1',       'LATENT_MRNA', 168163, 169474, '-'),
    # LATENT_BART (outside B95-8 deletion)
    ('BART_RPMS1',  'LATENT_BART', 151551, 153200, '+'),
    ('BART_BARF0',  'LATENT_BART', 159121, 160536, '+'),
    # POL3_EBER
    ('EBER-1',      'POL3_EBER',     6629,   6795, '+'),
    ('EBER-2',      'POL3_EBER',     6956,   7128, '+'),
    # LYTIC_IE
    ('BZLF1',       'LYTIC_IE',     89838,  90943, '-'),
    ('BRLF1',       'LYTIC_IE',     91078,  92894, '-'),
    # LYTIC_EARLY
    ('BHLF1',       'LYTIC_EARLY',  38287,  38380, '-'), # unique 5' part outside IR2
    ('BHRF1',       'LYTIC_EARLY',  42088,  42663, '+'),
    ('BFLF2',       'LYTIC_EARLY',  43691,  44647, '-'),
    ('BFLF1',       'LYTIC_EARLY',  44660,  46235, '-'),
    ('BFRF1A',      'LYTIC_EARLY',  46236,  46602, '+'),
    ('BFRF1',       'LYTIC_EARLY',  46603,  47519, '+'),
    ('BFRF2',       'LYTIC_EARLY',  47520,  49218, '+'),
    ('BORF2',       'LYTIC_EARLY',  64119,  66599, '+'),
    ('BaRF1',       'LYTIC_EARLY',  66612,  67520, '+'),
    ('BMRF1',       'LYTIC_EARLY',  67611,  68825, '+'),
    ('BMRF2',       'LYTIC_EARLY',  68830,  69903, '+'),
    ('BSLF2_BMLF1', 'LYTIC_EARLY',  70400,  71940, '-'),
    ('BSLF1',       'LYTIC_EARLY',  71969,  74593, '-'),
    ('BLLF3',       'LYTIC_EARLY',  75350,  76186, '-'),
    ('BRRF1',       'LYTIC_EARLY',  92895,  93826, '+'),
    ('BKRF3',       'LYTIC_EARLY',  98084,  98832, '+'),
    ('BBLF4',       'LYTIC_EARLY',  99542, 101915, '-'),
    ('BBLF2_BBLF3', 'LYTIC_EARLY', 105098, 106800, '-'),
    ('BGLF5',       'LYTIC_EARLY', 108687, 110036, '-'),
    ('BGLF4',       'LYTIC_EARLY', 110037, 111204, '-'),
    ('BGLF3.5',     'LYTIC_EARLY', 111205, 111652, '-'),
    ('BGLF3',       'LYTIC_EARLY', 111653, 112649, '-'),
    ('BDLF4',       'LYTIC_EARLY', 116087, 116711, '-'),
    ('BDLF3.5',     'LYTIC_EARLY', 116712, 116945, '-'),
    ('BcRF1',       'LYTIC_EARLY', 125179, 127415, '+'),
    ('BXLF1',       'LYTIC_EARLY', 130749, 132570, '-'),
    ('BVLF1',       'LYTIC_EARLY', 134840, 135627, '-'),
    ('BILF1',       'LYTIC_EARLY', 151703, 152641, '-'),
    ('BALF5',       'LYTIC_EARLY', 153241, 156288, '-'),
    ('BALF3',       'LYTIC_EARLY', 158865, 159120, '-'),
    ('BALF2',       'LYTIC_EARLY', 160926, 164312, '-'),
    ('BALF1',       'LYTIC_EARLY', 164397, 164945, '-'),
    ('BARF1',       'LYTIC_EARLY', 165046, 165711, '+'),
    ('BNLF2b',      'LYTIC_EARLY', 166540, 166836, '-'),
    ('BNLF2a',      'LYTIC_EARLY', 166846, 167028, '-'),
    # LYTIC_LATE (True Late / Structural)
    ('BNRF1',       'LYTIC_LATE',    1736,   5692, '+'),
    ('BCRF1',       'LYTIC_LATE',    9675,  10187, '+'),
    ('BFRF3',       'LYTIC_LATE',   49219,  49749, '+'),
    ('BPLF1',       'LYTIC_LATE',   49790,  59231, '-'),
    ('BOLF1',       'LYTIC_LATE',   59240,  62949, '-'),
    ('BORF1',       'LYTIC_LATE',   62950,  64044, '+'),
    ('BSRF1',       'LYTIC_LATE',   74636,  75292, '+'),
    ('BLRF1',       'LYTIC_LATE',   76259,  76567, '+'),
    ('BLRF2',       'LYTIC_LATE',   76637,  77125, '+'),
    ('BLLF1',       'LYTIC_LATE',   77142,  79865, '-'),
    ('BZLF2',       'LYTIC_LATE',   89157,  89828, '-'),
    ('BRRF2',       'LYTIC_LATE',   94014,  95627, '+'),
    ('BKRF2',       'LYTIC_LATE',   97670,  98064, '+'),
    ('BKRF4',       'LYTIC_LATE',   98846,  99499, '+'),
    ('BBRF1',       'LYTIC_LATE',  101916, 103659, '+'),
    ('BBRF2',       'LYTIC_LATE',  103660, 104496, '+'),
    ('BBRF3',       'LYTIC_LATE',  106849, 108066, '+'),
    ('BBLF1',       'LYTIC_LATE',  108459, 108640, '-'),
    ('BGRF1_BDRF1', 'LYTIC_LATE',  112650, 113574, '+'),
    ('BGLF2',       'LYTIC_LATE',  113575, 114562, '-'),
    ('BGLF1',       'LYTIC_LATE',  114563, 116055, '-'),
    ('BDLF3',       'LYTIC_LATE',  118074, 118778, '-'),
    ('BDLF2',       'LYTIC_LATE',  118839, 120101, '-'),
    ('BDLF1',       'LYTIC_LATE',  120112, 121017, '-'),
    ('BcLF1',       'LYTIC_LATE',  121033, 125178, '-'),
    ('BTRF1',       'LYTIC_LATE',  127416, 128626, '+'),
    ('BXLF2',       'LYTIC_LATE',  128627, 130747, '-'),
    ('BXRF1',       'LYTIC_LATE',  132571, 133126, '+'),
    ('BVRF1',       'LYTIC_LATE',  133127, 134808, '+'),
    ('BVRF2_BdRF1', 'LYTIC_LATE',  135638, 137455, '+'),
    ('BILF2',       'LYTIC_LATE',  137490, 138236, '-'),
    ('BALF4',       'LYTIC_LATE',  156291, 158850, '-'),
    # B958_DEL (Wild-type EBV RNA detector)
    ('LF3',         'B958_DEL',    140570, 143344, '-'),
    ('LF2',         'B958_DEL',    149035, 150284, '-'),
    ('LF1',         'B958_DEL',    150285, 151400, '-'),
]

gene_counts = {g[0]: 0 for g in GENES}
gene_len    = {g[0]: (g[3] - g[2] + 1) for g in GENES}
gene_class  = {g[0]: g[1] for g in GENES}
class_len   = {}
for g, k, s, e, _ in GENES:
    class_len[k] = class_len.get(k, 0) + (e - s + 1)

cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

total_ebv_mapped = 0
host_decoy_filtered = 0
bcrf1_host_filtered = 0
barf1_host_filtered = 0
ir1_mapq0_reads = 0
bzlf1_spliced_reads = 0
ebna_spliced_reads = 0
lmp_spliced_reads = 0

for line in open(sam_path):
    if line.startswith('@'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 11: continue
    qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
    if flag & 2308: continue # skip unmapped, secondary, supplementary

    # Parse CIGAR for aligned reference span, matched bases, and splice junctions ('N')
    ref_span = 0
    aligned_bp = 0
    splices = []
    cur_pos = pos
    for length_s, op in cigar_re.findall(cigar):
        l = int(length_s)
        if op in ('M', '=', 'X'):
            ref_span += l; aligned_bp += l; cur_pos += l
        elif op == 'D':
            ref_span += l; cur_pos += l
        elif op == 'N':
            if l >= 40:
                splices.append((cur_pos, cur_pos + l - 1))
            ref_span += l; cur_pos += l

    if aligned_bp < 35: continue
    mid = pos + ref_span // 2

    # Competitive host filter: did this read match a GENCODE v45 human transcript >= 90% of its EBV aligned_bp?
    qkey = qname + ("/1" if (flag & 64) else ("/2" if (flag & 128) else ""))
    h_match = max(host_best.get(qkey, 0), host_best.get(qname, 0))
    if h_match >= 45 and h_match >= 0.90 * aligned_bp:
        host_decoy_filtered += 1
        if 9600 <= mid <= 10300: bcrf1_host_filtered += 1
        if 165000 <= mid <= 165800: barf1_host_filtered += 1
        continue

    # Track IR1 multi-mapping reads (12,301-35,037) before MAPQ filter
    if 12301 <= mid <= 35037:
        ir1_mapq0_reads += 1
        # Check for canonical EBNA W1-W2 / Cp / Wp splice junctions
        if any(11300 <= j_s <= 36500 for j_s, _ in splices):
            ebna_spliced_reads += 1
        continue

    if mapq < 10: continue
    total_ebv_mapped += 1

    # Splice junction classification outside IR1
    for j_s, j_e in splices:
        if 89900 <= j_s <= 90850 and (j_e - j_s) <= 250:
            bzlf1_spliced_reads += 1
        elif 11000 <= j_s <= 97500:
            ebna_spliced_reads += 1
        elif j_s >= 166000 or j_e <= 1700:
            lmp_spliced_reads += 1

    # Assign read by midpoint to curated non-overlapping gene interval
    for g, k, s, e, _ in GENES:
        if s <= mid <= e:
            gene_counts[g] += 1
            break

# Write per-gene table
with open(out_genes, 'w') as gh:
    gh.write('run\tdonor\tgene\tkinetic_class\tlength_bp\treads\trpm\trpkm\n')
    for g, k, s, e, _ in GENES:
        cnt = gene_counts[g]
        glen = e - s + 1
        rpm  = cnt * 1e6 / num_reads if num_reads > 0 else 0.0
        rpkm = rpm * 1000.0 / glen if glen > 0 else 0.0
        gh.write(f'{run_id}\t{donor}\t{g}\t{k}\t{glen}\t{cnt}\t{rpm:.4f}\t{rpkm:.4f}\n')

# Aggregate by kinetic class
class_counts = {}
for g, cnt in gene_counts.items():
    k = gene_class[g]
    class_counts[k] = class_counts.get(k, 0) + cnt

def get_rpkm(k):
    cnt = class_counts.get(k, 0)
    glen = class_len.get(k, 1)
    return (cnt * 1e6 / num_reads) * (1000.0 / glen) if num_reads > 0 else 0.0

def get_rpm(k):
    return class_counts.get(k, 0) * 1e6 / num_reads if num_reads > 0 else 0.0

lat_rpkm   = get_rpkm('LATENT_MRNA')
bart_rpkm  = get_rpkm('LATENT_BART')
eber_rpkm  = get_rpkm('POL3_EBER')
ie_rpkm    = get_rpkm('LYTIC_IE')
early_rpkm = get_rpkm('LYTIC_EARLY')
late_rpkm  = get_rpkm('LYTIC_LATE')
del_rpkm   = get_rpkm('B958_DEL')

# Primary lytic score: log2((late_rpkm + 0.05) / (lat_rpkm + 0.05)) and exact ratio
lytic_late_lat_ratio = (late_rpkm / lat_rpkm) if lat_rpkm > 0 else 0.0
lytic_late_lat_log2  = math.log2((late_rpkm + 0.05) / (lat_rpkm + 0.05))
lytic_ie_lat_log2    = math.log2((ie_rpkm + 0.05) / (lat_rpkm + 0.05))
lytic_early_lat_log2 = math.log2((early_rpkm + 0.05) / (lat_rpkm + 0.05))

cols = [
    'run', 'donor', 'population', 'continentalGroup', 'batch', 'RIN', 'total_reads',
    'ebv_unique_reads', 'ir1_mapq0_reads', 'host_decoy_filtered', 'bcrf1_host_filtered', 'barf1_host_filtered',
    'latent_mrna_reads', 'latent_bart_reads', 'pol3_eber_reads',
    'lytic_ie_reads', 'lytic_early_reads', 'lytic_late_reads', 'b958_del_reads',
    'bzlf1_spliced_reads', 'ebna_spliced_reads', 'lmp_spliced_reads',
    'latent_mrna_rpm', 'lytic_ie_rpm', 'lytic_early_rpm', 'lytic_late_rpm', 'b958_del_rpm',
    'latent_mrna_rpkm', 'latent_bart_rpkm', 'pol3_eber_rpkm',
    'lytic_ie_rpkm', 'lytic_early_rpkm', 'lytic_late_rpkm', 'b958_del_rpkm',
    'lytic_late_over_latent_ratio', 'lytic_late_over_latent_log2',
    'lytic_ie_over_latent_log2', 'lytic_early_over_latent_log2'
]
vals = [
    run_id, donor, pop, cont, batch, rin, f'{int(num_reads)}',
    f'{total_ebv_mapped}', f'{ir1_mapq0_reads}', f'{host_decoy_filtered}', f'{bcrf1_host_filtered}', f'{barf1_host_filtered}',
    f'{class_counts.get("LATENT_MRNA",0)}', f'{class_counts.get("LATENT_BART",0)}', f'{class_counts.get("POL3_EBER",0)}',
    f'{class_counts.get("LYTIC_IE",0)}', f'{class_counts.get("LYTIC_EARLY",0)}', f'{class_counts.get("LYTIC_LATE",0)}', f'{class_counts.get("B958_DEL",0)}',
    f'{bzlf1_spliced_reads}', f'{ebna_spliced_reads}', f'{lmp_spliced_reads}',
    f'{get_rpm("LATENT_MRNA"):.4f}', f'{get_rpm("LYTIC_IE"):.4f}', f'{get_rpm("LYTIC_EARLY"):.4f}', f'{get_rpm("LYTIC_LATE"):.4f}', f'{get_rpm("B958_DEL"):.4f}',
    f'{lat_rpkm:.4f}', f'{bart_rpkm:.4f}', f'{eber_rpkm:.4f}',
    f'{ie_rpkm:.4f}', f'{early_rpkm:.4f}', f'{late_rpkm:.4f}', f'{del_rpkm:.4f}',
    f'{lytic_late_lat_ratio:.5f}', f'{lytic_late_lat_log2:.4f}',
    f'{lytic_ie_lat_log2:.4f}', f'{lytic_early_lat_log2:.4f}'
]
with open(out_tsv, 'w') as oh:
    oh.write('\t'.join(cols) + '\n')
    oh.write('\t'.join(vals) + '\n')
