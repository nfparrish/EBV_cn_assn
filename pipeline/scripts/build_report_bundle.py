import csv
import json
import math
import pathlib
import sys


def to_float(value):
    try:
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return None
        return out
    except Exception:
        return None


def load_rows(path):
    with pathlib.Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def top_rows(rows, source, max_rows=25):
    selected = []
    for row in rows:
        p = to_float(row.get('p_value'))
        fdr = to_float(row.get('p_value_fdr'))
        effect = to_float(row.get('effect_size'))
        feature_id = row.get('variant_id') or row.get('feature_id') or row.get('gene_id') or ''
        selected.append({
            'source': source,
            'pathway': row.get('pathway', ''),
            'layer': row.get('layer', ''),
            'feature_id': feature_id,
            'gene_id': row.get('gene_id', ''),
            'effect_size': effect,
            'p_value': p,
            'p_value_fdr': fdr,
            'n_samples': row.get('n_samples', ''),
        })

    selected.sort(key=lambda x: (
        x['p_value'] if x['p_value'] is not None else 1.0,
        x['p_value_fdr'] if x['p_value_fdr'] is not None else 1.0,
        -abs(x['effect_size']) if x['effect_size'] is not None else 0.0,
    ))
    return selected[:max_rows]


def main():
    assoc_diag_path = pathlib.Path(sys.argv[1])
    assoc_qc_summary_path = pathlib.Path(sys.argv[2])
    assoc_qc_hits_path = pathlib.Path(sys.argv[3])
    rna_qc_summary_path = pathlib.Path(sys.argv[4])
    rna_qc_hits_path = pathlib.Path(sys.argv[5])
    rna_variant_presence_path = pathlib.Path(sys.argv[6])
    rna_ingest_report_path = pathlib.Path(sys.argv[7])

    assoc_diag = json.loads(assoc_diag_path.read_text())
    assoc_qc_summary = json.loads(assoc_qc_summary_path.read_text())
    rna_qc_summary = json.loads(rna_qc_summary_path.read_text())
    rna_variant_presence = json.loads(rna_variant_presence_path.read_text())
    rna_ingest = json.loads(rna_ingest_report_path.read_text())

    assoc_hits = load_rows(assoc_qc_hits_path)
    rna_hits = load_rows(rna_qc_hits_path)

    merged = top_rows(assoc_hits, 'wgs_assoc', max_rows=20) + top_rows(rna_hits, 'rna_assoc', max_rows=20)
    merged.sort(key=lambda x: (
        x['p_value'] if x['p_value'] is not None else 1.0,
        x['p_value_fdr'] if x['p_value_fdr'] is not None else 1.0,
    ))

    with pathlib.Path('report_bundle_top_signals.tsv').open('w') as handle:
        fields = ['source', 'pathway', 'layer', 'feature_id', 'gene_id', 'effect_size', 'p_value', 'p_value_fdr', 'n_samples']
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in merged[:40]:
            writer.writerow({
                'source': row['source'],
                'pathway': row['pathway'],
                'layer': row['layer'],
                'feature_id': row['feature_id'],
                'gene_id': row['gene_id'],
                'effect_size': '' if row['effect_size'] is None else f"{row['effect_size']:.8g}",
                'p_value': '' if row['p_value'] is None else f"{row['p_value']:.8g}",
                'p_value_fdr': '' if row['p_value_fdr'] is None else f"{row['p_value_fdr']:.8g}",
                'n_samples': row['n_samples'],
            })

    metrics = {
        'wgs_assoc': {
            'ebv_model_coverage_rate': assoc_diag.get('ebv', {}).get('model_coverage_rate'),
            'host_model_coverage_rate': assoc_diag.get('host', {}).get('model_coverage_rate'),
            'overall_fdr_lt_0_10_total': assoc_diag.get('overall', {}).get('fdr_lt_0_10_total'),
            'qc_top_hit_count': assoc_qc_summary.get('top_hit_count'),
        },
        'rna_assoc': {
            'ebv_gene_count': rna_ingest.get('ebv_gene_count'),
            'model_health': rna_qc_summary.get('model_health', {}),
            'variant_presence_flags': rna_variant_presence.get('overall_flags', {}),
        },
        'bundle': {
            'merged_top_signal_count': min(len(merged), 40),
        },
    }
    pathlib.Path('report_bundle_metrics.json').write_text(json.dumps(metrics, indent=2))

    summary_lines = [
        '# EBV Multi-Omics Reporting Bundle',
        '',
        '## WGS Association Highlights',
        f"- EBV model coverage rate: {metrics['wgs_assoc']['ebv_model_coverage_rate']}",
        f"- Host model coverage rate: {metrics['wgs_assoc']['host_model_coverage_rate']}",
        f"- Overall FDR<0.10 signals: {metrics['wgs_assoc']['overall_fdr_lt_0_10_total']}",
        f"- WGS QC top-hit count: {metrics['wgs_assoc']['qc_top_hit_count']}",
        '',
        '## RNA Integration Highlights',
        f"- EBV genes detected in RNA ingest: {metrics['rna_assoc']['ebv_gene_count']}",
        f"- RNA model health: {json.dumps(metrics['rna_assoc']['model_health'])}",
        f"- Variant presence flags: {json.dumps(metrics['rna_assoc']['variant_presence_flags'])}",
        '',
        '## Bundle Artifacts',
        f"- merged_top_signal_count: {metrics['bundle']['merged_top_signal_count']}",
        '- Top signals: report_bundle_top_signals.tsv',
        '- Machine-readable metrics: report_bundle_metrics.json',
    ]
    pathlib.Path('report_bundle_summary.md').write_text('\n'.join(summary_lines) + '\n')


if __name__ == '__main__':
    main()
