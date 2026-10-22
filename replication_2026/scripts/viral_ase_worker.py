#!/usr/bin/env python3
import csv, os, re, sys
from collections import Counter

run_id, sam_path, out_ase, out_spl = sys.argv[1:5]
BASE = '/work/nfp8/EBV'

# Load reference chrEBV sequence (1-based indexing via seq[pos-1])
ref_seq = ''.join(l.strip() for l in open('/hpc/group/parrishlab/refs/ebv/chrEBV.fa') if not l.startswith('>')).upper()

# Load all candidate sites called in this donor (or across the cohort)
# We will record base/indel pileup at all positions in variants/all_calls.tsv that have >=5 carriers
site_set = set()
for line in open(f'{BASE}/variants/pon_sites.tsv'):
    if line.startswith('pos'): continue
    f = line.rstrip('\n').split('\t')
    if int(f[3]) >= 5:
        site_set.add((int(f[0]), f[1], f[2]))

cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
# For each position p in site_set, count (ref_matches, alt_matches, total_cov)
pos_to_sites = {}
for p, r, a in site_set:
    pos_to_sites.setdefault(p, []).append((r, a))

counts = {s: [0, 0] for s in site_set} # [ref_cnt, alt_cnt]

# Fine-grained splice counters
spl_counts = Counter()

for line in open(sam_path):
    if line.startswith('@'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 11: continue
    flag, pos, mapq, cigar, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
    if flag & 2308: continue

    # Walk CIGAR to map reference positions to read offsets and extract splice junctions ('N')
    r_pos = pos
    q_pos = 0
    ops = [(int(l), op) for l, op in cigar_re.findall(cigar)]
    for idx, (l, op) in enumerate(ops):
        if op in ('M', '=', 'X'):
            # Check if any target site falls in [r_pos, r_pos + l - 1]
            if mapq >= 10:
                for p in range(r_pos, r_pos + l):
                    if p in pos_to_sites:
                        q_idx = q_pos + (p - r_pos)
                        # Check if next op is I or D for indel sites at p
                        next_op = ops[idx+1] if (p == r_pos + l - 1 and idx + 1 < len(ops)) else None
                        for ref_a, alt_a in pos_to_sites[p]:
                            if len(ref_a) == 1 and len(alt_a) == 1:
                                if next_op and next_op[1] in ('I', 'D'): continue
                                b = seq[q_idx]
                                if b == ref_a: counts[(p, ref_a, alt_a)][0] += 1
                                elif b == alt_a: counts[(p, ref_a, alt_a)][1] += 1
                            elif len(alt_a) > len(ref_a): # Insertion after p
                                ins_len = len(alt_a) - len(ref_a)
                                if next_op and next_op[1] == 'I' and next_op[0] == ins_len:
                                    counts[(p, ref_a, alt_a)][1] += 1
                                elif not next_op or next_op[1] not in ('I', 'D'):
                                    counts[(p, ref_a, alt_a)][0] += 1
                            elif len(ref_a) > len(alt_a): # Deletion after p
                                del_len = len(ref_a) - len(alt_a)
                                if next_op and next_op[1] == 'D' and next_op[0] == del_len:
                                    counts[(p, ref_a, alt_a)][1] += 1
                                elif not next_op or next_op[1] not in ('I', 'D'):
                                    counts[(p, ref_a, alt_a)][0] += 1
            r_pos += l; q_pos += l
        elif op == 'I':
            q_pos += l
        elif op == 'D':
            r_pos += l
        elif op == 'S':
            q_pos += l
        elif op == 'N':
            if l >= 40:
                j_s, j_e = r_pos, r_pos + l - 1
                # Classify canonical EBV splice junctions
                if 11300 <= j_s <= 12000:
                    spl_counts['Cp_EBNA'] += 1
                elif 12001 <= j_s <= 35100 and 12001 <= j_e <= 35355:
                    spl_counts['Wp_W1_W2_IR1'] += 1
                elif 12001 <= j_s <= 35355 and 35356 <= j_e <= 36300:
                    spl_counts['W2_to_Y1Y2_EBNA2'] += 1
                elif 12001 <= j_s <= 36500 and 41900 <= j_e <= 42200:
                    spl_counts['Latent_W2_Y_to_BHRF1'] += 1
                elif 62300 <= j_s <= 62600 and 95500 <= j_e <= 95800:
                    spl_counts['Qp_to_EBNA1'] += 1
                elif 36000 <= j_s <= 68000 and 95500 <= j_e <= 95800:
                    spl_counts['CpWp_U_to_EBNA1'] += 1
                elif 90500 <= j_s <= 90700 and 90500 <= j_e <= 90700:
                    spl_counts['BZLF1_ex1_ex2'] += 1
                elif 90200 <= j_s <= 90350 and 90200 <= j_e <= 90350:
                    spl_counts['BZLF1_ex2_ex3'] += 1
                elif 166400 <= j_s <= 166600 or (j_s <= 1700 and j_e <= 1700):
                    spl_counts['LMP2A_2B_exons'] += 1
                elif 168100 <= j_s <= 169500:
                    spl_counts['LMP1_introns'] += 1
                elif 151500 <= j_s <= 160600:
                    spl_counts['BART_splices'] += 1
            r_pos += l

with open(out_ase, 'w') as oh:
    oh.write('run\tpos\tref\talt\trna_ref\trna_alt\n')
    for (p, r, a), (rc, ac) in sorted(counts.items()):
        if rc + ac > 0:
            oh.write(f'{run_id}\t{p}\t{r}\t{a}\t{rc}\t{ac}\n')

with open(out_spl, 'w') as sh:
    sh.write('run\tsplice_class\tcount\n')
    for k in ['Cp_EBNA', 'Wp_W1_W2_IR1', 'W2_to_Y1Y2_EBNA2', 'Latent_W2_Y_to_BHRF1', 'Qp_to_EBNA1', 'CpWp_U_to_EBNA1', 'BZLF1_ex1_ex2', 'BZLF1_ex2_ex3', 'LMP2A_2B_exons', 'LMP1_introns', 'BART_splices']:
        sh.write(f'{run_id}\t{k}\t{spl_counts[k]}\n')
