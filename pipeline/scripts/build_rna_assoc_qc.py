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


def load_tsv_rows(path):
    rows = []
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def summarize_assoc(rows):
    total = len(rows)
    finite_p = 0
    nominal = 0
    fdr_hits = 0
    modeled = 0

    for row in rows:
        p = to_float(row.get('p_value'))
        fdr = to_float(row.get('p_value_fdr'))
        effect = to_float(row.get('effect_size'))

        if p is not None:
            finite_p += 1
            if p < 0.05:
                nominal += 1
        if fdr is not None and fdr < 0.10:
            fdr_hits += 1
        if effect is not None:
            modeled += 1

    finite_rate = (finite_p / total) if total else 0.0
    modeled_rate = (modeled / total) if total else 0.0

    return {
        'test_count': total,
        'finite_p_count': finite_p,
        'finite_p_rate': round(finite_rate, 6),
        'modeled_count': modeled,
        'modeled_rate': round(modeled_rate, 6),
        'nominal_hit_count': nominal,
        'fdr_hit_count': fdr_hits,
    }


def top_hits(rows, label, max_hits=10):
    ranked = []
    for row in rows:
        p = to_float(row.get('p_value'))
        effect = to_float(row.get('effect_size'))
        if p is None:
            continue
        ranked.append({
            'layer': label,
            'feature_id': row.get('variant_id', row.get('gene_id', '')),
            'gene_id': row.get('gene_id', ''),
            'effect_size': effect,
            'p_value': p,
            'p_value_fdr': to_float(row.get('p_value_fdr')),
            'n_samples': row.get('n_samples', ''),
        })

    ranked.sort(key=lambda x: (x['p_value'], -abs(x['effect_size']) if x['effect_size'] is not None else 0.0))
    return ranked[:max_hits]


def main():
    cn_log2_results = pathlib.Path(sys.argv[1])
    cn_log2_report = pathlib.Path(sys.argv[2])
    cn_quant_results = pathlib.Path(sys.argv[3])
    cn_quant_report = pathlib.Path(sys.argv[4])
    var_log2_results = pathlib.Path(sys.argv[5])
    var_log2_report = pathlib.Path(sys.argv[6])
    var_quant_results = pathlib.Path(sys.argv[7])
    var_quant_report = pathlib.Path(sys.argv[8])

    cn_log2_rows = load_tsv_rows(cn_log2_results)
    cn_quant_rows = load_tsv_rows(cn_quant_results)
    var_log2_rows = load_tsv_rows(var_log2_results)
    var_quant_rows = load_tsv_rows(var_quant_results)

    reports = {
        'cn_log2': json.loads(cn_log2_report.read_text()),
        'cn_quantile': json.loads(cn_quant_report.read_text()),
        'variant_log2': json.loads(var_log2_report.read_text()),
        'variant_quantile': json.loads(var_quant_report.read_text()),
    }

    summary = {
        'coverage': {
            'sample_overlap': {
                'cn_log2': reports['cn_log2'].get('sample_overlap_count', 0),
                'cn_quantile': reports['cn_quantile'].get('sample_overlap_count', 0),
                'variant_log2': reports['variant_log2'].get('sample_overlap_count', 0),
                'variant_quantile': reports['variant_quantile'].get('sample_overlap_count', 0),
            },
            'gene_counts': {
                'cn_log2': reports['cn_log2'].get('gene_count', 0),
                'cn_quantile': reports['cn_quantile'].get('gene_count', 0),
                'variant_log2': reports['variant_log2'].get('gene_count', 0),
                'variant_quantile': reports['variant_quantile'].get('gene_count', 0),
            },
            'variant_count_modeled': {
                'variant_log2': reports['variant_log2'].get('variant_count_modeled', 0),
                'variant_quantile': reports['variant_quantile'].get('variant_count_modeled', 0),
            },
        },
        'significance': {
            'cn_log2': summarize_assoc(cn_log2_rows),
            'cn_quantile': summarize_assoc(cn_quant_rows),
            'variant_log2': summarize_assoc(var_log2_rows),
            'variant_quantile': summarize_assoc(var_quant_rows),
        },
    }

    all_stats = [
        summary['significance']['cn_log2'],
        summary['significance']['cn_quantile'],
        summary['significance']['variant_log2'],
        summary['significance']['variant_quantile'],
    ]

    summary['model_health'] = {
        'all_layers_have_finite_p': all(stat['finite_p_count'] > 0 for stat in all_stats),
        'low_coverage_any_layer': any(stat['finite_p_rate'] < 0.5 for stat in all_stats),
        'all_layers_nominally_null': all(stat['nominal_hit_count'] == 0 for stat in all_stats),
    }

    hits = []
    hits.extend(top_hits(cn_log2_rows, 'cn_log2'))
    hits.extend(top_hits(cn_quant_rows, 'cn_quantile'))
    hits.extend(top_hits(var_log2_rows, 'variant_log2'))
    hits.extend(top_hits(var_quant_rows, 'variant_quantile'))
    hits.sort(key=lambda x: (x['p_value'], -abs(x['effect_size']) if x['effect_size'] is not None else 0.0))

    with pathlib.Path('rna_assoc_qc_top_hits.tsv').open('w') as handle:
        header = ['layer', 'feature_id', 'gene_id', 'effect_size', 'p_value', 'p_value_fdr', 'n_samples']
        handle.write('\t'.join(header) + '\n')
        for hit in hits[:40]:
            handle.write('\t'.join([
                hit['layer'],
                str(hit['feature_id']),
                str(hit['gene_id']),
                str(hit['effect_size']),
                str(hit['p_value']),
                str(hit['p_value_fdr']),
                str(hit['n_samples']),
            ]) + '\n')

    pathlib.Path('rna_assoc_qc_summary.json').write_text(json.dumps(summary, indent=2))

    report_lines = [
        '# RNA Association QC Report',
        '',
        '## Coverage',
        f"- CN log2 sample overlap: {summary['coverage']['sample_overlap']['cn_log2']}",
        f"- CN quantile sample overlap: {summary['coverage']['sample_overlap']['cn_quantile']}",
        f"- Variant log2 sample overlap: {summary['coverage']['sample_overlap']['variant_log2']}",
        f"- Variant quantile sample overlap: {summary['coverage']['sample_overlap']['variant_quantile']}",
        '',
        '## Significance',
        f"- CN log2 nominal hits (p<0.05): {summary['significance']['cn_log2']['nominal_hit_count']}",
        f"- CN quantile nominal hits (p<0.05): {summary['significance']['cn_quantile']['nominal_hit_count']}",
        f"- Variant log2 nominal hits (p<0.05): {summary['significance']['variant_log2']['nominal_hit_count']}",
        f"- Variant quantile nominal hits (p<0.05): {summary['significance']['variant_quantile']['nominal_hit_count']}",
        '',
        '## Model Health',
        f"- all_layers_have_finite_p: {summary['model_health']['all_layers_have_finite_p']}",
        f"- low_coverage_any_layer: {summary['model_health']['low_coverage_any_layer']}",
        f"- all_layers_nominally_null: {summary['model_health']['all_layers_nominally_null']}",
        '',
        'Top hits are available in rna_assoc_qc_top_hits.tsv.',
    ]

    pathlib.Path('rna_assoc_qc_report.md').write_text('\n'.join(report_lines) + '\n')


if __name__ == '__main__':
    main()
