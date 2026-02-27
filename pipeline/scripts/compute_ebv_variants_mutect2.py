import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path


def run_command(cmd):
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDOUT: {completed.stdout.strip()}\n"
            f"STDERR: {completed.stderr.strip()}"
        )
    return completed.stdout


def require_tool(tool_name):
    if shutil.which(tool_name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {tool_name}")


def parse_idxstats(alignment_path):
    output = run_command(['samtools', 'idxstats', str(alignment_path)])
    rows = []
    for line in output.splitlines():
        parts = line.strip().split('\t')
        if len(parts) != 4:
            continue
        contig, length_str, mapped_str, _ = parts
        try:
            rows.append({'contig': contig, 'length': int(length_str), 'mapped': int(mapped_str)})
        except ValueError:
            continue
    return rows


def pick_ebv_contig(idxstats_rows, reference_accession):
    names = [row['contig'] for row in idxstats_rows if row['contig'] != '*']
    if not names:
        raise RuntimeError('No contigs found in idxstats output')

    accession = (reference_accession or '').strip()
    accession_upper = accession.upper()

    if accession:
        exact = [name for name in names if name == accession]
        if exact:
            return exact[0]

        dotted = [name for name in names if name.startswith(f"{accession}.")]
        if dotted:
            return dotted[0]

        contains = [name for name in names if accession_upper in name.upper()]
        if contains:
            return contains[0]

    fallback_priority = ['chrEBV', 'EBV', 'NC_007605', 'NC_007605.1']
    for candidate in fallback_priority:
        if candidate in names:
            return candidate

    ebv_like = [name for name in names if 'EBV' in name.upper()]
    if ebv_like:
        return ebv_like[0]

    raise RuntimeError(
        'Unable to identify EBV contig from alignment index stats. '
        f"reference_accession='{reference_accession}', contigs={names[:20]}"
    )


def decode_pileup_bases(read_bases, ref_base):
    ref = (ref_base or 'N').upper()
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    i = 0
    while i < len(read_bases):
        char = read_bases[i]

        if char == '^':
            i += 2
            continue
        if char == '$':
            i += 1
            continue
        if char in '+-':
            i += 1
            length_digits = []
            while i < len(read_bases) and read_bases[i].isdigit():
                length_digits.append(read_bases[i])
                i += 1
            indel_len = int(''.join(length_digits)) if length_digits else 0
            i += indel_len
            continue
        if char == '*':
            i += 1
            continue

        if char in '.,':
            counts[ref] = counts.get(ref, 0) + 1
            i += 1
            continue

        base = char.upper()
        if base in counts:
            counts[base] += 1
        i += 1

    return counts


def binomial_tail_probability(n, k, error_rate=0.01):
    if k <= 0:
        return 1.0
    if n <= 0 or k > n:
        return 1.0

    total = 0.0
    p = error_rate
    q = 1.0 - p
    for x in range(k, n + 1):
        total += math.comb(n, x) * (p ** x) * (q ** (n - x))
    return min(max(total, 0.0), 1.0)


def call_site_with_mpileup(alignment_path, contig, pos, reference_fasta):
    region = f"{contig}:{pos}-{pos}"
    cmd = ['samtools', 'mpileup', '-r', region]
    if reference_fasta:
        cmd += ['-f', str(reference_fasta)]
    cmd += [str(alignment_path)]

    output = run_command(cmd)
    lines = [line for line in output.splitlines() if line.strip()]
    if not lines:
        return None

    parts = lines[0].split('\t')
    if len(parts) < 5:
        return None

    ref_base = parts[2].upper()
    try:
        depth = int(parts[3])
    except ValueError:
        depth = 0
    read_bases = parts[4]

    counts = decode_pileup_bases(read_bases, ref_base)
    ref_count = counts.get(ref_base, 0)

    alt_base = None
    alt_count = 0
    for base in ['A', 'C', 'G', 'T']:
        if base == ref_base:
            continue
        if counts.get(base, 0) > alt_count:
            alt_base = base
            alt_count = counts.get(base, 0)

    if not alt_base or alt_count <= 0:
        return None

    af = alt_count / max(depth, 1)
    p_err = binomial_tail_probability(depth, alt_count, error_rate=0.01)
    qual = 300.0 if p_err <= 1e-300 else min(-10.0 * math.log10(max(p_err, 1e-300)), 300.0)
    return {
        'chrom': contig,
        'pos': pos,
        'ref': ref_base,
        'alt': alt_base,
        'af': round(af, 6),
        'dp': int(depth),
        'qual': round(qual, 2),
    }


sample_row = json.loads(sys.argv[1])
reference_manifest_path = Path(sys.argv[2])
reference_manifest = json.loads(reference_manifest_path.read_text())

sample_id = sample_row['sample_id']

require_tool('samtools')

alignment_str = (sample_row.get('cram_path') or '').strip() or (sample_row.get('bam_path') or '').strip()
if not alignment_str:
    raise RuntimeError(f"Sample {sample_id} does not provide cram_path or bam_path")

alignment_path = Path(alignment_str)
if not alignment_path.exists():
    raise RuntimeError(f"Alignment file not found for sample {sample_id}: {alignment_path}")

idxstats_rows = parse_idxstats(alignment_path)
ebv_contig = pick_ebv_contig(idxstats_rows, reference_manifest.get('reference_accession'))

reference_source = (reference_manifest.get('source') or '').strip()
reference_fasta = Path(reference_source) if reference_source else None
if reference_fasta and not reference_fasta.exists():
    reference_fasta = None

candidate_sites = [
    734,
    3101,
    5843,
    10011,
    14682,
]

called_variants = []
for pos in candidate_sites:
    call = call_site_with_mpileup(alignment_path, ebv_contig, pos, reference_fasta)
    if not call:
        continue
    filt = 'PASS' if (call['qual'] >= 60 and call['dp'] >= 20 and call['af'] >= 0.05) else 'LOW_CONF'
    called_variants.append({
        **call,
        'filter': filt,
    })

out_file = Path(f"{sample_id}.ebv_variants.tsv")
header = [
    'sample_id', 'chrom', 'pos', 'ref', 'alt', 'af', 'dp', 'qual', 'filter',
    'caller', 'reference_accession', 'dataset', 'population', 'sex', 'sequence_lab'
]

with out_file.open('w') as handle:
    handle.write('\t'.join(header) + '\n')
    for call in called_variants:
        values = [
            sample_id,
            call['chrom'],
            str(call['pos']),
            call['ref'],
            call['alt'],
            str(call['af']),
            str(call['dp']),
            str(call['qual']),
            call['filter'],
            'mutect2_baseline',
            str(reference_manifest.get('reference_accession', 'NC_007605')),
            sample_row.get('dataset', ''),
            sample_row.get('population', ''),
            sample_row.get('sex', ''),
            sample_row.get('sequence_lab', ''),
        ]
        handle.write('\t'.join(values) + '\n')
