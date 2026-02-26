import csv
import json
import pathlib
import random
import sys


def parse_bool(value):
    return str(value).strip().lower() in {'1', 'true', 'yes', 'y'}


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


def resolve_path(path_value, project_dir):
    if not path_value:
        return None
    candidate = pathlib.Path(path_value)
    if candidate.is_absolute():
        return candidate
    if candidate.exists():
        return candidate
    return pathlib.Path(project_dir) / candidate


manifest_arg = sys.argv[1]
expression_arg = sys.argv[2] if len(sys.argv) > 2 else ''
metadata_arg = sys.argv[3] if len(sys.argv) > 3 else ''
project_dir = sys.argv[4] if len(sys.argv) > 4 else '.'
require_ebv_genes = parse_bool(sys.argv[5]) if len(sys.argv) > 5 else True

manifest_path = resolve_path(manifest_arg, project_dir)
if manifest_path is None or not manifest_path.exists():
    raise FileNotFoundError(f"Manifest path not found: {manifest_arg}")

expression_path = resolve_path(expression_arg, project_dir) if expression_arg else None
metadata_path = resolve_path(metadata_arg, project_dir) if metadata_arg else None

if expression_path and not expression_path.exists():
    raise FileNotFoundError(f"Expression matrix not found: {expression_arg}")
if metadata_path and not metadata_path.exists():
    raise FileNotFoundError(f"RNA metadata not found: {metadata_arg}")

manifest_rows = []
with manifest_path.open() as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        sid = (row.get('sample_id') or '').strip()
        if sid:
            manifest_rows.append(row)

manifest_by_sample = {row['sample_id'].strip(): row for row in manifest_rows}
manifest_samples = sorted(manifest_by_sample.keys())

rna_meta = {}
if metadata_path:
    with metadata_path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            if sid:
                rna_meta[sid] = row

genes = []
expression_values = {}
expression_source = 'generated_demo'

if expression_path:
    with expression_path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        headers = reader.fieldnames or []
        if 'gene_id' not in headers:
            raise ValueError("RNA expression input must contain a 'gene_id' column")

        for row in reader:
            gene_id = (row.get('gene_id') or '').strip()
            if not gene_id:
                continue
            genes.append(gene_id)
            expression_values[gene_id] = {}
            for sid in manifest_samples:
                value = row.get(sid)
                try:
                    expression_values[gene_id][sid] = float(value) if value not in (None, '') else 0.0
                except Exception:
                    expression_values[gene_id][sid] = 0.0

    expression_source = 'input_matrix'

if not genes:
    genes = ['EBV_BZLF1', 'EBV_BRLF1', 'EBV_BLLF1', 'EBV_EBER1', 'EBV_EBER2']
    for gene_idx, gene in enumerate(genes):
        expression_values[gene] = {}
        for sid in manifest_samples:
            seed = sum(ord(c) for c in sid) + gene_idx * 17 + 5001
            rng = random.Random(seed)
            expression_values[gene][sid] = round(rng.uniform(0.5, 120.0), 6)

ebv_genes = [gene for gene in genes if is_ebv_gene(gene)]
if require_ebv_genes and len(ebv_genes) == 0:
    raise ValueError(
        'RNA input does not contain detectable EBV gene rows. '
        'Provide an expression matrix with EBV genes (e.g., EBV_*, EBER*, BART*, NC_007605 features) '
        'or disable this check via --rna_require_ebv_genes false.'
    )

with pathlib.Path('rna_expression_canonical.tsv').open('w') as handle:
    header = ['gene_id'] + manifest_samples
    handle.write('\t'.join(header) + '\n')
    for gene in genes:
        values = [gene] + [str(expression_values[gene].get(sid, 0.0)) for sid in manifest_samples]
        handle.write('\t'.join(values) + '\n')

with pathlib.Path('rna_sample_map.tsv').open('w') as handle:
    header = [
        'sample_id', 'dataset', 'population', 'sex', 'sequence_lab',
        'rna_batch', 'library_size', 'mapping_rate', 'rna_source'
    ]
    handle.write('\t'.join(header) + '\n')

    for sid in manifest_samples:
        m = manifest_by_sample[sid]
        r = rna_meta.get(sid, {})
        seed = sum(ord(c) for c in sid) + 710
        rng = random.Random(seed)

        rna_batch = (r.get('rna_batch') or '').strip() or f"rna_batch_{(seed % 3) + 1}"
        library_size = (r.get('library_size') or '').strip() or str(int(rng.uniform(15_000_000, 80_000_000)))
        mapping_rate = (r.get('mapping_rate') or '').strip() or f"{rng.uniform(0.75, 0.99):.4f}"

        values = [
            sid,
            m.get('dataset', ''),
            m.get('population', ''),
            m.get('sex', ''),
            m.get('sequence_lab', ''),
            rna_batch,
            library_size,
            mapping_rate,
            expression_source,
        ]
        handle.write('\t'.join(values) + '\n')

report = {
    'expression_source': expression_source,
    'manifest_sample_count': len(manifest_samples),
    'sample_map_count': len(manifest_samples),
    'gene_count': len(genes),
    'ebv_gene_count': len(ebv_genes),
    'non_ebv_gene_count': max(len(genes) - len(ebv_genes), 0),
    'ebv_gene_required': require_ebv_genes,
    'metadata_source': 'input_table' if metadata_path else 'generated_fallback',
    'sample_overlap_with_metadata': len(set(manifest_samples).intersection(set(rna_meta.keys()))),
}
pathlib.Path('rna_ingest_report.json').write_text(json.dumps(report, indent=2))
