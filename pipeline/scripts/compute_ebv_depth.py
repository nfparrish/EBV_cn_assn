import json
import shutil
import subprocess
import sys
from pathlib import Path


def run_command(cmd, extra_env=None):
    import os
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    completed = subprocess.run(cmd, capture_output=True, text=True, env=env)
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


def ref_args(cram_reference, tool='view'):
    """Return reference flag args appropriate for the samtools tool.
    flagstat does not accept a reference flag; use CRAM_REFERENCE env var instead."""
    ref = (cram_reference or '').strip()
    if not ref:
        return []
    if tool == 'coverage':
        return ['--reference', ref]
    if tool == 'view':
        return ['-T', ref]
    # flagstat: no CLI flag supported; caller must set CRAM_REFERENCE env var
    return []


def parse_header_contigs(alignment_path, cram_reference):
    """Get contig names from the CRAM/BAM header using samtools view -H.
    This is fast (reads only the header, not alignment records)."""
    cmd = ['samtools', 'view', '-H'] + ref_args(cram_reference, 'view') + [str(alignment_path)]
    output = run_command(cmd)
    contigs = []
    for line in output.splitlines():
        if line.startswith('@SQ'):
            parts = dict(f.split(':', 1) for f in line.split('\t')[1:] if ':' in f)
            name = parts.get('SN')
            if name:
                contigs.append(name)
    return contigs


def pick_ebv_contig(contig_names, reference_accession):
    names = [n for n in contig_names if n != '*']
    if not names:
        raise RuntimeError('No contigs found in alignment header')

    accession = (reference_accession or '').strip()
    accession_upper = accession.upper()

    if accession:
        exact = [n for n in names if n == accession]
        if exact:
            return exact[0]
        dotted = [n for n in names if n.startswith(f"{accession}.")]
        if dotted:
            return dotted[0]
        contains = [n for n in names if accession_upper in n.upper()]
        if contains:
            return contains[0]

    fallback_priority = ['chrEBV', 'EBV', 'NC_007605', 'NC_007605.1']
    for candidate in fallback_priority:
        if candidate in names:
            return candidate

    ebv_like = [n for n in names if 'EBV' in n.upper()]
    if ebv_like:
        return ebv_like[0]

    raise RuntimeError(
        'Unable to identify EBV contig from alignment header. '
        f"reference_accession='{reference_accession}', contigs={names[:20]}"
    )


def pick_host_depth_contig(contig_names):
    """Choose a single large autosomal contig as the host depth baseline.
    Using one chromosome is sufficient for normalization and avoids 22 sequential
    samtools coverage calls. chr1 (or '1') is preferred: longest autosome, stable depth."""
    available = set(contig_names)
    for candidate in ['chr1', '1']:
        if candidate in available:
            return candidate
    # Fallback: first autosomal contig found
    chr_names = [f'chr{i}' for i in range(1, 23)] + [str(i) for i in range(1, 23)]
    for c in chr_names:
        if c in available:
            return c
    raise RuntimeError('No autosomal contig found for host depth baseline')


def get_coverage_stats(alignment_path, cram_reference, contigs):
    """Run samtools coverage for the given contig list.
    Returns dict: contig -> {'meandepth': float, 'numreads': int}"""
    region_args = []
    for c in contigs:
        region_args += ['-r', c]
    cmd = (['samtools', 'coverage'] +
           ref_args(cram_reference, 'coverage') +
           region_args + [str(alignment_path)])
    output = run_command(cmd)
    stats = {}
    for line in output.splitlines():
        if not line or line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        rname = parts[0]
        try:
            numreads = int(parts[3])
            meandepth = float(parts[6])
        except ValueError:
            continue
        stats[rname] = {'meandepth': meandepth, 'numreads': numreads}
    return stats


def weighted_mean_depth(stats, contigs):
    """Compute length-weighted mean depth across contigs from coverage stats."""
    # samtools coverage reports meandepth directly; we average by numreads weight
    # Use simple unweighted mean of meandepth values (sufficient approximation).
    depths = [stats[c]['meandepth'] for c in contigs if c in stats]
    if not depths:
        return 0.0
    return sum(depths) / len(depths)


def total_mapped_reads(alignment_path, cram_reference, threads=4):
    """Get total mapped reads via samtools flagstat.
    flagstat does not support a --reference flag; relies on CRAM_REFERENCE env var."""
    cmd = ['samtools', 'flagstat', f'-@{threads}', str(alignment_path)]
    output = run_command(cmd)
    for line in output.splitlines():
        if 'mapped (' in line and 'primary' not in line:
            return int(line.split()[0])
    # fallback: parse "in total" line
    for line in output.splitlines():
        if 'in total' in line:
            return int(line.split()[0])
    return 0


row = json.loads(sys.argv[1])
reference_manifest_path = Path(sys.argv[2])
reference = json.loads(reference_manifest_path.read_text())
cram_reference = sys.argv[3].strip() if len(sys.argv) > 3 else ''
sample_id = row['sample_id']

require_tool('samtools')

alignment_str = (row.get('cram_path') or '').strip() or (row.get('bam_path') or '').strip()
if not alignment_str:
    raise RuntimeError(f"Sample {sample_id} does not provide cram_path or bam_path")

alignment_path = Path(alignment_str)
if not alignment_path.exists():
    raise RuntimeError(f"Alignment file not found for sample {sample_id}: {alignment_path}")

contig_names = parse_header_contigs(alignment_path, cram_reference)
ebv_contig = pick_ebv_contig(contig_names, reference.get('reference_accession'))
host_contig = pick_host_depth_contig(contig_names)

# Get coverage stats: EBV and host baseline contig in separate calls
ebv_stats = get_coverage_stats(alignment_path, cram_reference, [ebv_contig])
host_stats = get_coverage_stats(alignment_path, cram_reference, [host_contig])

ebv_depth = round(ebv_stats.get(ebv_contig, {}).get('meandepth', 0.0), 6)
ebv_reads = ebv_stats.get(ebv_contig, {}).get('numreads', 0)
host_depth = round(host_stats.get(host_contig, {}).get('meandepth', 0.0), 6)
total_reads_value = total_mapped_reads(alignment_path, cram_reference)

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
    'host_depth_contig': host_contig,
    'metrics': {
        'host_autosomal_depth': host_depth,
        'ebv_depth': ebv_depth,
        'total_reads': total_reads_value,
        'ebv_mapped_reads': ebv_reads,
    }
}

Path(f"{sample_id}.ebv_depth.json").write_text(json.dumps(payload, indent=2))
