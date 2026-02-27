import json
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
        contig, length_str, mapped_str, unmapped_str = parts
        try:
            rows.append({
                'contig': contig,
                'length': int(length_str),
                'mapped': int(mapped_str),
                'unmapped': int(unmapped_str),
            })
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


def autosomal_contigs(idxstats_rows):
    available = {row['contig'] for row in idxstats_rows}
    chr_names = [f'chr{i}' for i in range(1, 23)]
    plain_names = [str(i) for i in range(1, 23)]
    if any(name in available for name in chr_names):
        return [name for name in chr_names if name in available]
    return [name for name in plain_names if name in available]


def parse_coverage_row(line):
    parts = line.strip().split('\t')
    if len(parts) < 7:
        return None
    try:
        start = int(parts[1])
        end = int(parts[2])
        meandepth = float(parts[6])
    except ValueError:
        return None
    length = max(end - start + 1, 0)
    return {'length': length, 'meandepth': meandepth}


def mean_depth_for_contigs(alignment_path, contigs):
    total_bases = 0
    total_depth = 0.0
    for contig in contigs:
        output = run_command(['samtools', 'coverage', '-r', contig, str(alignment_path)])
        for line in output.splitlines():
            if not line or line.startswith('#'):
                continue
            parsed = parse_coverage_row(line)
            if not parsed:
                continue
            total_bases += parsed['length']
            total_depth += parsed['meandepth'] * parsed['length']

    if total_bases == 0:
        return 0.0
    return total_depth / total_bases


def mapped_reads_for_contig(idxstats_rows, contig):
    for row in idxstats_rows:
        if row['contig'] == contig:
            return max(row['mapped'], 0)
    return 0


def total_reads(alignment_path):
    output = run_command(['samtools', 'view', '-c', str(alignment_path)])
    return int(output.strip() or 0)


row = json.loads(sys.argv[1])
reference_manifest_path = Path(sys.argv[2])
reference = json.loads(reference_manifest_path.read_text())
sample_id = row['sample_id']

require_tool('samtools')

alignment_str = (row.get('cram_path') or '').strip() or (row.get('bam_path') or '').strip()
if not alignment_str:
    raise RuntimeError(f"Sample {sample_id} does not provide cram_path or bam_path")

alignment_path = Path(alignment_str)
if not alignment_path.exists():
    raise RuntimeError(f"Alignment file not found for sample {sample_id}: {alignment_path}")

idxstats_rows = parse_idxstats(alignment_path)
ebv_contig = pick_ebv_contig(idxstats_rows, reference.get('reference_accession'))
host_contigs = autosomal_contigs(idxstats_rows)
if not host_contigs:
    raise RuntimeError(
        f"No autosomal contigs found for sample {sample_id}. "
        "Expected chr1..chr22 or 1..22 in alignment header/index."
    )

host_depth = round(mean_depth_for_contigs(alignment_path, host_contigs), 6)
ebv_depth = round(mean_depth_for_contigs(alignment_path, [ebv_contig]), 6)
total_reads_value = total_reads(alignment_path)
ebv_reads = mapped_reads_for_contig(idxstats_rows, ebv_contig)

payload = {
    'sample_id': sample_id,
    'dataset': row.get('dataset', ''),
    'population': row.get('population', ''),
    'sex': row.get('sex', ''),
    'sequence_lab': row.get('sequence_lab', ''),
    'reference_accession': reference.get('reference_accession'),
    'reference_source_type': reference.get('source_type'),
    'alignment_path': str(alignment_path),
    'alignment_type': 'cram' if (row.get('cram_path') or '').strip() else 'bam',
    'ebv_contig': ebv_contig,
    'metrics': {
        'host_autosomal_depth': host_depth,
        'ebv_depth': ebv_depth,
        'total_reads': total_reads_value,
        'ebv_mapped_reads': ebv_reads,
    }
}

Path(f"{sample_id}.ebv_depth.json").write_text(json.dumps(payload, indent=2))
