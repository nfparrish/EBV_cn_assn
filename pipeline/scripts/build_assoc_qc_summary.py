import csv
import json
import pathlib
import sys


def to_float(value):
    try:
        return float(value)
    except Exception:
        return None


def load_rows(path):
    rows = []
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def extract_hits(rows, pathway):
    hits = []
    for row in rows:
        p = to_float(row.get('p_value'))
        fdr = to_float(row.get('p_value_fdr'))
        effect = to_float(row.get('effect_size'))
        if p is None:
            continue
        if p < 0.05 or (fdr is not None and fdr < 0.10):
            hits.append({
                'variant_id': row.get('variant_id', ''),
                'pathway': pathway,
                'effect_size': effect,
                'p_value': p,
                'p_value_fdr': fdr,
                'n_samples': row.get('n_samples', ''),
            })
    hits.sort(key=lambda x: (x['p_value'] if x['p_value'] is not None else 1.0, abs(x['effect_size']) if x['effect_size'] is not None else 0.0), reverse=False)
    return hits


diagnostics = json.loads(pathlib.Path(sys.argv[1]).read_text())
ebv_results = load_rows(sys.argv[2])
host_results = load_rows(sys.argv[3])

hits = extract_hits(ebv_results, 'ebv') + extract_hits(host_results, 'host')
hits.sort(key=lambda x: (x['p_value'] if x['p_value'] is not None else 1.0))

with pathlib.Path('association_qc_top_hits.tsv').open('w') as handle:
    header = ['pathway', 'variant_id', 'effect_size', 'p_value', 'p_value_fdr', 'n_samples']
    handle.write('\t'.join(header) + '\n')
    for hit in hits[:20]:
        handle.write('\t'.join([
            hit['pathway'],
            hit['variant_id'],
            str(hit['effect_size']),
            str(hit['p_value']),
            str(hit['p_value_fdr']),
            str(hit['n_samples']),
        ]) + '\n')

qc_summary = {
    'ebv': diagnostics.get('ebv', {}),
    'host': diagnostics.get('host', {}),
    'overall': diagnostics.get('overall', {}),
    'top_hit_count': min(len(hits), 20),
    'high_level_flags': {
        'low_model_coverage_any_pathway': any(
            (diagnostics.get(k, {}).get('model_coverage_rate', 0.0) < 0.5)
            for k in ['ebv', 'host']
        ),
        'no_nominal_hits': len(hits) == 0,
    },
}

pathlib.Path('association_qc_summary.json').write_text(json.dumps(qc_summary, indent=2))

report_lines = [
    '# Association QC Report',
    '',
    '## Coverage',
    f"- EBV modeled variants: {diagnostics.get('ebv', {}).get('variant_count_modeled')} / {diagnostics.get('ebv', {}).get('variant_count_total')}",
    f"- Host modeled variants: {diagnostics.get('host', {}).get('variant_count_modeled')} / {diagnostics.get('host', {}).get('variant_count_total')}",
    '',
    '## Signal summary',
    f"- EBV p<0.05: {diagnostics.get('ebv', {}).get('p_lt_0_05_count')}",
    f"- Host p<0.05: {diagnostics.get('host', {}).get('p_lt_0_05_count')}",
    f"- Overall FDR<0.10: {diagnostics.get('overall', {}).get('fdr_lt_0_10_total')}",
    '',
    '## Flags',
    f"- low_model_coverage_any_pathway: {qc_summary['high_level_flags']['low_model_coverage_any_pathway']}",
    f"- no_nominal_hits: {qc_summary['high_level_flags']['no_nominal_hits']}",
    '',
    'Top hits are available in association_qc_top_hits.tsv.',
]
pathlib.Path('association_qc_report.md').write_text('\n'.join(report_lines) + '\n')
