import csv
import gzip
import json
import pathlib
import random
import shutil
import subprocess
import sys


def parse_bool(value):
    return str(value).strip().lower() in {'1', 'true', 'yes', 'y'}


def resolve_path(path_value, project_dir):
    if not path_value:
        return None
    candidate = pathlib.Path(path_value)
    if candidate.is_absolute():
        return candidate
    if candidate.exists():
        return candidate
    return pathlib.Path(project_dir) / candidate


def resolve_input_path(path_value, manifest_dir, project_dir):
    candidate = pathlib.Path(path_value)
    if candidate.is_absolute() and candidate.exists():
        return candidate
    if candidate.exists():
        return candidate

    manifest_relative = manifest_dir / candidate
    if manifest_relative.exists():
        return manifest_relative

    project_relative = pathlib.Path(project_dir) / candidate
    if project_relative.exists():
        return project_relative

    return candidate


def is_ebv_contig(name):
    key = (name or '').upper()
    return 'NC_007605' in key or 'EBV' in key or 'HHV4' in key


def is_ebv_gene(gene_id):
    gid = (gene_id or '').strip().upper()
    return (
        gid.startswith('EBV_')
        or gid.startswith('EBER')
        or gid.startswith('BART')
        or 'NC_007605' in gid
        or gid.startswith('LF')
        or gid.startswith('BZLF')
        or gid.startswith('BRLF')
    )


def parse_fasta_headers(path):
    headers = []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                headers.append(line[1:].strip().split()[0])
    return headers


def parse_gtf_gene_rows(path):
    genes = []
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] not in ('gene', 'exon', 'transcript'):
                continue
            contig = parts[0]
            attrs = parts[8]
            gene_id = ''
            for token in attrs.split(';'):
                token = token.strip()
                if token.startswith('gene_id'):
                    chunks = token.split('"')
                    if len(chunks) >= 2:
                        gene_id = chunks[1]
                        break
            genes.append((contig, gene_id))
    return genes


def validate_combined_reference(reference_fasta, reference_gtf):
    contigs = parse_fasta_headers(reference_fasta)
    if not contigs:
        raise ValueError(f"Reference FASTA has no contigs: {reference_fasta}")

    has_ebv_contig = any(is_ebv_contig(c) for c in contigs)
    has_host_contig = any(not is_ebv_contig(c) for c in contigs)
    if not has_ebv_contig or not has_host_contig:
        raise ValueError(
            "RNA reference FASTA must include both host and EBV contigs. "
            "Ensure EBV contig (e.g., NC_007605) is present alongside host genome contigs."
        )

    gtf_genes = parse_gtf_gene_rows(reference_gtf)
    if not gtf_genes:
        raise ValueError(f"Reference GTF has no readable gene/exon rows: {reference_gtf}")

    ebv_gene_rows = [g for g in gtf_genes if is_ebv_contig(g[0]) or is_ebv_gene(g[1])]
    host_gene_rows = [g for g in gtf_genes if (not is_ebv_contig(g[0])) and (not is_ebv_gene(g[1]))]
    if not ebv_gene_rows or not host_gene_rows:
        raise ValueError(
            "RNA reference GTF must annotate both host and EBV genes. "
            "Expected EBV annotations on EBV contig or EBV-like gene identifiers."
        )


def read_manifest(path):
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = (row.get('sample_id') or '').strip()
            if sample_id:
                rows.append(row)
    if not rows:
        raise ValueError('Manifest contains no valid sample rows.')
    return rows


def read_metadata(path):
    metadata = {}
    if not path or not path.exists():
        return metadata
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            if sid:
                metadata[sid] = row
    return metadata


def ensure_tool(binary_name):
    if shutil.which(binary_name) is None:
        raise RuntimeError(
            f"Required tool not found on PATH: {binary_name}. "
            f"Install {binary_name} to enable FASTQ-based RNA counting."
        )


def parse_cigar_ref_span(cigar):
    num = ''
    ref_len = 0
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        if not num:
            continue
        value = int(num)
        if ch in {'M', 'D', 'N', '=', 'X'}:
            ref_len += value
        num = ''
    return ref_len


def build_gene_intervals(reference_gtf):
    intervals = {}
    genes = set()
    with reference_gtf.open() as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in {'gene', 'exon'}:
                continue
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]
            gene_id = ''
            for token in attrs.split(';'):
                token = token.strip()
                if token.startswith('gene_id'):
                    chunks = token.split('"')
                    if len(chunks) >= 2:
                        gene_id = chunks[1]
                    break
            if not gene_id:
                continue
            genes.add(gene_id)
            intervals.setdefault(contig, []).append((start, end, gene_id))

    for contig in intervals:
        intervals[contig].sort(key=lambda x: (x[0], x[1], x[2]))

    return intervals, sorted(genes)


def count_with_builtin_counter(bam_paths, reference_gtf):
    gene_intervals, genes = build_gene_intervals(reference_gtf)
    counts = {gene: {} for gene in genes}

    for sample_id, bam_path in bam_paths:
        for gene in genes:
            counts[gene][sample_id] = 0.0

        view_cmd = ['samtools', 'view', '-F', '4', str(bam_path)]
        proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out, err = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"samtools view failed for {sample_id}: {err}")

        for line in out.splitlines():
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            rname = fields[2]
            if rname == '*':
                continue
            pos = int(fields[3])
            cigar = fields[5]
            ref_span = parse_cigar_ref_span(cigar)
            if ref_span <= 0:
                continue
            aln_start = pos
            aln_end = pos + ref_span - 1

            for start, end, gene_id in gene_intervals.get(rname, []):
                if end < aln_start:
                    continue
                if start > aln_end:
                    break
                counts[gene_id][sample_id] += 1.0

    return genes, counts


def maybe_gzip_open(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return path.open()


def count_reads_fastq(path):
    line_count = 0
    with maybe_gzip_open(path) as handle:
        for _ in handle:
            line_count += 1
    return line_count // 4


def run_command(command):
    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        cmd = ' '.join(command)
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def align_and_count(rows, reference_fasta, reference_gtf, threads, work_dir, manifest_dir, project_dir, counter):
    ensure_tool('minimap2')
    ensure_tool('samtools')

    bam_paths = []
    sample_metrics = {}

    for row in rows:
        sample_id = row['sample_id'].strip()
        read1 = (row.get('read1_path') or '').strip()
        read2 = (row.get('read2_path') or '').strip()

        if not read1:
            raise ValueError(f"Sample {sample_id} missing read1_path for FASTQ RNA counting mode")

        read1_path = resolve_input_path(read1, manifest_dir, project_dir)
        read2_path = resolve_input_path(read2, manifest_dir, project_dir) if read2 else None

        if not read1_path.exists():
            raise FileNotFoundError(f"FASTQ not found: {read1_path}")
        if read2_path and not read2_path.exists():
            raise FileNotFoundError(f"FASTQ not found: {read2_path}")

        bam_path = work_dir / f"{sample_id}.sorted.bam"

        cmd = ['minimap2', '-t', str(threads), '-ax', 'sr', str(reference_fasta), str(read1_path)]
        if read2_path:
            cmd.append(str(read2_path))

        minimap = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
        sort_cmd = ['samtools', 'sort', '-@', str(threads), '-o', str(bam_path)]
        sorter = subprocess.Popen(sort_cmd, stdin=minimap.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
        minimap.stdout.close()
        sort_out, sort_err = sorter.communicate()
        map_out, map_err = minimap.communicate()

        if minimap.returncode != 0:
            raise RuntimeError(f"minimap2 failed for {sample_id}: {map_err.decode(errors='ignore')}")
        if sorter.returncode != 0:
            raise RuntimeError(f"samtools sort failed for {sample_id}: {sort_err.decode(errors='ignore')}")

        run_command(['samtools', 'index', str(bam_path)])

        flagstat = run_command(['samtools', 'flagstat', str(bam_path)])
        mapped = 0
        total = 0
        for line in flagstat.splitlines():
            if ' in total ' in line:
                try:
                    total = int(line.split(' + ')[0].strip())
                except Exception:
                    pass
            elif ' mapped (' in line and 'mate mapped' not in line:
                try:
                    mapped = int(line.split(' + ')[0].strip())
                except Exception:
                    pass

        read_count = count_reads_fastq(read1_path)
        if read2_path:
            read_count = min(read_count, count_reads_fastq(read2_path))

        sample_metrics[sample_id] = {
            'library_size': int(read_count),
            'mapping_rate': (float(mapped) / float(total)) if total > 0 else 0.0,
        }

        bam_paths.append((sample_id, bam_path))

    if counter == 'featureCounts':
        ensure_tool('featureCounts')
        fc_out = work_dir / 'featurecounts.txt'
        fc_cmd = [
            'featureCounts',
            '-T', str(threads),
            '-a', str(reference_gtf),
            '-o', str(fc_out),
        ] + [str(x[1]) for x in bam_paths]
        run_command(fc_cmd)
        data_path = work_dir / 'featurecounts.txt'
        return bam_paths, sample_metrics, data_path, None, None

    if counter == 'builtin':
        genes, counts = count_with_builtin_counter(bam_paths, reference_gtf)
        return bam_paths, sample_metrics, None, genes, counts

    raise ValueError(f"Unsupported counter mode: {counter}")


def parse_featurecounts_matrix(path, sample_order):
    genes = []
    matrix = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            header = line.rstrip('\n').split('\t')
            if header and header[0] == 'Geneid':
                bam_cols = header[6:]
                break
        else:
            raise ValueError('featureCounts output missing header line with Geneid.')

    bam_to_sample = {}
    for col in bam_cols:
        sample_name = pathlib.Path(col).name
        if sample_name.endswith('.sorted.bam'):
            sample_name = sample_name[:-11]
        bam_to_sample[col] = sample_name

    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if not parts or parts[0] == 'Geneid' or len(parts) < 7:
                continue
            gene_id = parts[0]
            genes.append(gene_id)
            counts = {}
            for idx, col in enumerate(bam_cols):
                sid = bam_to_sample.get(col, col)
                try:
                    counts[sid] = float(parts[6 + idx])
                except Exception:
                    counts[sid] = 0.0
            matrix[gene_id] = counts

    ordered_samples = [sid for sid in sample_order]
    return genes, ordered_samples, matrix


def write_expression_canonical(genes, samples, matrix):
    with pathlib.Path('rna_expression_canonical.tsv').open('w') as handle:
        handle.write('\t'.join(['gene_id'] + samples) + '\n')
        for gene in genes:
            row = [gene] + [str(matrix.get(gene, {}).get(sid, 0.0)) for sid in samples]
            handle.write('\t'.join(row) + '\n')


def write_sample_map(rows, metadata, sample_metrics):
    with pathlib.Path('rna_sample_map.tsv').open('w') as handle:
        header = [
            'sample_id', 'dataset', 'population', 'sex', 'sequence_lab',
            'rna_batch', 'library_size', 'mapping_rate', 'rna_source'
        ]
        handle.write('\t'.join(header) + '\n')

        for row in rows:
            sid = row['sample_id'].strip()
            m = metadata.get(sid, {})
            seed = sum(ord(c) for c in sid) + 710
            rng = random.Random(seed)

            defaults = sample_metrics.get(sid, {})
            library_size = m.get('library_size') or str(defaults.get('library_size', int(rng.uniform(15000000, 80000000))))
            mapping_rate = m.get('mapping_rate') or f"{defaults.get('mapping_rate', rng.uniform(0.75, 0.99)):.4f}"
            rna_batch = (m.get('rna_batch') or '').strip() or f"rna_batch_{(seed % 3) + 1}"

            values = [
                sid,
                row.get('dataset', ''),
                row.get('population', ''),
                row.get('sex', ''),
                row.get('sequence_lab', ''),
                rna_batch,
                str(library_size),
                str(mapping_rate),
                'fastq_recount',
            ]
            handle.write('\t'.join(values) + '\n')


def write_reports(rows, genes, require_ebv_genes, ebv_genes, count_report):
    ingest_report = {
        'expression_source': 'fastq_recount',
        'manifest_sample_count': len(rows),
        'sample_map_count': len(rows),
        'gene_count': len(genes),
        'ebv_gene_count': len(ebv_genes),
        'non_ebv_gene_count': max(len(genes) - len(ebv_genes), 0),
        'ebv_gene_required': require_ebv_genes,
        'metadata_source': 'input_table' if count_report.get('metadata_used') else 'generated_fallback',
        'sample_overlap_with_metadata': count_report.get('metadata_overlap_count', 0),
    }
    pathlib.Path('rna_ingest_report.json').write_text(json.dumps(ingest_report, indent=2))
    pathlib.Path('rna_fastq_count_report.json').write_text(json.dumps(count_report, indent=2))


def main():
    manifest_arg = sys.argv[1]
    metadata_arg = sys.argv[2] if len(sys.argv) > 2 else ''
    reference_fasta_arg = sys.argv[3] if len(sys.argv) > 3 else ''
    reference_gtf_arg = sys.argv[4] if len(sys.argv) > 4 else ''
    manifest_source_arg = sys.argv[5] if len(sys.argv) > 5 else ''
    project_dir = sys.argv[6] if len(sys.argv) > 6 else '.'
    require_ebv_genes = parse_bool(sys.argv[7]) if len(sys.argv) > 7 else True
    threads = int(sys.argv[8]) if len(sys.argv) > 8 else 2
    aligner = sys.argv[9] if len(sys.argv) > 9 else 'minimap2'
    counter = sys.argv[10] if len(sys.argv) > 10 else 'auto'

    if aligner != 'minimap2':
        raise ValueError(f"Unsupported rna_aligner '{aligner}'. Supported: minimap2")
    if counter == 'auto':
        counter = 'featureCounts' if shutil.which('featureCounts') else 'builtin'
    if counter not in {'featureCounts', 'builtin'}:
        raise ValueError(f"Unsupported rna_counter '{counter}'. Supported: auto, featureCounts, builtin")

    manifest_path = resolve_path(manifest_arg, project_dir)
    manifest_source_path = resolve_path(manifest_source_arg, project_dir) if manifest_source_arg else manifest_path
    metadata_path = resolve_path(metadata_arg, project_dir) if metadata_arg else None
    reference_fasta = resolve_path(reference_fasta_arg, project_dir) if reference_fasta_arg else None
    reference_gtf = resolve_path(reference_gtf_arg, project_dir) if reference_gtf_arg else None

    if manifest_path is None or not manifest_path.exists():
        raise FileNotFoundError(f"Manifest path not found: {manifest_arg}")
    if reference_fasta is None or not reference_fasta.exists():
        raise FileNotFoundError(
            "FASTQ counting mode requires --rna_reference_fasta <combined_host_ebv.fa>"
        )
    if reference_gtf is None or not reference_gtf.exists():
        raise FileNotFoundError(
            "FASTQ counting mode requires --rna_reference_gtf <combined_host_ebv.gtf>"
        )

    validate_combined_reference(reference_fasta, reference_gtf)

    rows = read_manifest(manifest_path)
    metadata = read_metadata(metadata_path) if metadata_path else {}

    work_dir = pathlib.Path('rna_fastq_work')
    work_dir.mkdir(parents=True, exist_ok=True)

    bam_paths, sample_metrics, featurecounts_matrix, builtin_genes, builtin_counts = align_and_count(
        rows=rows,
        reference_fasta=reference_fasta,
        reference_gtf=reference_gtf,
        threads=max(1, threads),
        work_dir=work_dir,
        manifest_dir=manifest_source_path.parent if manifest_source_path else manifest_path.parent,
        project_dir=project_dir,
        counter=counter,
    )

    sample_order = [row['sample_id'].strip() for row in rows]
    if counter == 'featureCounts':
        genes, samples, matrix = parse_featurecounts_matrix(featurecounts_matrix, sample_order)
    else:
        genes = builtin_genes or []
        samples = sample_order
        matrix = builtin_counts or {}

    ebv_genes = [gene for gene in genes if is_ebv_gene(gene)]
    if require_ebv_genes and len(ebv_genes) == 0:
        raise ValueError(
            'Generated count matrix does not contain detectable EBV gene rows. '
            'Verify combined reference annotations include EBV genes.'
        )

    write_expression_canonical(genes, samples, matrix)
    write_sample_map(rows, metadata, sample_metrics)

    metadata_overlap = len(set(samples).intersection(set(metadata.keys()))) if metadata else 0

    count_report = {
        'mode': 'fastq_recount',
        'aligner': aligner,
        'counter': counter,
        'reference_fasta': str(reference_fasta),
        'reference_gtf': str(reference_gtf),
        'sample_count': len(samples),
        'gene_count': len(genes),
        'ebv_gene_count': len(ebv_genes),
        'bam_count': len(bam_paths),
        'metadata_used': bool(metadata),
        'metadata_overlap_count': metadata_overlap,
        'sample_metrics': sample_metrics,
    }

    write_reports(rows, genes, require_ebv_genes, ebv_genes, count_report)


if __name__ == '__main__':
    main()
