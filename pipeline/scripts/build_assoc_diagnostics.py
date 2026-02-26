import csv
import json
import math
import pathlib
import statistics
import sys


def to_float(value):
    try:
        v = float(value)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    except Exception:
        return None


def load_results(path):
    rows = []
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def summarize(pathway, rows):
    total = len(rows)
    finite_effect = [to_float(r.get('effect_size')) for r in rows]
    finite_effect = [x for x in finite_effect if x is not None]

    finite_se = [to_float(r.get('std_error')) for r in rows]
    finite_se = [x for x in finite_se if x is not None]

    finite_p = [to_float(r.get('p_value')) for r in rows]
    finite_p = [x for x in finite_p if x is not None]

    finite_fdr = [to_float(r.get('p_value_fdr')) for r in rows]
    finite_fdr = [x for x in finite_fdr if x is not None]

    z_values = [to_float(r.get('z_score')) for r in rows]
    z_values = [x for x in z_values if x is not None]

    low_p = sum(1 for x in finite_p if x < 0.05)
    low_fdr = sum(1 for x in finite_fdr if x < 0.10)

    return {
        'pathway': pathway,
        'variant_count_total': total,
        'variant_count_modeled': len(finite_p),
        'variant_count_unmodeled': total - len(finite_p),
        'effect_abs_median': statistics.median([abs(x) for x in finite_effect]) if finite_effect else None,
        'std_error_median': statistics.median(finite_se) if finite_se else None,
        'z_abs_median': statistics.median([abs(x) for x in z_values]) if z_values else None,
        'p_lt_0_05_count': low_p,
        'fdr_lt_0_10_count': low_fdr,
        'model_coverage_rate': (len(finite_p) / total) if total else 0.0,
        'convergence_proxy': 'modeled_if_finite_effect_and_pvalue',
        'residual_proxy': 'z_score_distribution_summary',
    }


ebv_path = pathlib.Path(sys.argv[1])
host_path = pathlib.Path(sys.argv[2])

ebv_rows = load_results(ebv_path)
host_rows = load_results(host_path)

ebv_summary = summarize('ebv', ebv_rows)
host_summary = summarize('host', host_rows)

combined = {
    'ebv': ebv_summary,
    'host': host_summary,
    'overall': {
        'variant_count_total': ebv_summary['variant_count_total'] + host_summary['variant_count_total'],
        'variant_count_modeled': ebv_summary['variant_count_modeled'] + host_summary['variant_count_modeled'],
        'variant_count_unmodeled': ebv_summary['variant_count_unmodeled'] + host_summary['variant_count_unmodeled'],
        'fdr_lt_0_10_total': ebv_summary['fdr_lt_0_10_count'] + host_summary['fdr_lt_0_10_count'],
    }
}

pathlib.Path('association_diagnostics_summary.json').write_text(json.dumps(combined, indent=2))

with pathlib.Path('association_diagnostics.tsv').open('w') as handle:
    header = [
        'pathway', 'variant_count_total', 'variant_count_modeled', 'variant_count_unmodeled',
        'effect_abs_median', 'std_error_median', 'z_abs_median',
        'p_lt_0_05_count', 'fdr_lt_0_10_count', 'model_coverage_rate',
    ]
    handle.write('\t'.join(header) + '\n')

    for item in [ebv_summary, host_summary]:
        values = [
            item['pathway'],
            str(item['variant_count_total']),
            str(item['variant_count_modeled']),
            str(item['variant_count_unmodeled']),
            str(item['effect_abs_median']),
            str(item['std_error_median']),
            str(item['z_abs_median']),
            str(item['p_lt_0_05_count']),
            str(item['fdr_lt_0_10_count']),
            str(round(item['model_coverage_rate'], 6)),
        ]
        handle.write('\t'.join(values) + '\n')
