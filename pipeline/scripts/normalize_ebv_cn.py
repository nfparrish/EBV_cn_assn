import glob
import json
import pathlib
import statistics

rows = [json.loads(pathlib.Path(p).read_text()) for p in sorted(glob.glob('*.ebv_cn_raw.json'))]

depth_values = [r['methods']['depth_ratio_cn'] for r in rows] if rows else []
read_values = [r['methods']['read_count_ratio_cn'] for r in rows] if rows else []

depth_mean = statistics.mean(depth_values) if depth_values else 1.0
read_mean = statistics.mean(read_values) if read_values else 1.0

def safe_div(v, denom):
    return 0.0 if denom == 0 else v / denom

raw_header = [
    'sample_id','dataset','population','sex','sequence_lab','reference_accession',
    'depth_ratio_cn','read_count_ratio_cn','host_autosomal_depth','ebv_depth','total_reads','ebv_mapped_reads'
]
with pathlib.Path('ebv_cn_raw_cohort.tsv').open('w') as handle:
    handle.write('\t'.join(raw_header) + '\n')
    for row in rows:
        vals = [
            row['sample_id'], row.get('dataset',''), row.get('population',''), row.get('sex',''), row.get('sequence_lab',''),
            str(row.get('reference_accession','')),
            str(row['methods']['depth_ratio_cn']),
            str(row['methods']['read_count_ratio_cn']),
            str(row['support_metrics']['host_autosomal_depth']),
            str(row['support_metrics']['ebv_depth']),
            str(row['support_metrics']['total_reads']),
            str(row['support_metrics']['ebv_mapped_reads'])
        ]
        handle.write('\t'.join(vals) + '\n')

norm_header = [
    'sample_id','dataset','population','sex','sequence_lab','reference_accession',
    'depth_ratio_cn','read_count_ratio_cn','depth_ratio_cn_norm','read_count_ratio_cn_norm'
]
with pathlib.Path('ebv_cn_normalized.tsv').open('w') as handle:
    handle.write('\t'.join(norm_header) + '\n')
    for row in rows:
        d = row['methods']['depth_ratio_cn']
        r = row['methods']['read_count_ratio_cn']
        vals = [
            row['sample_id'], row.get('dataset',''), row.get('population',''), row.get('sex',''), row.get('sequence_lab',''),
            str(row.get('reference_accession','')),
            str(d),
            str(r),
            str(round(safe_div(d, depth_mean), 6)),
            str(round(safe_div(r, read_mean), 6)),
        ]
        handle.write('\t'.join(vals) + '\n')

summary = {
    'sample_count': len(rows),
    'methods': {
        'depth_ratio_cn_mean': round(depth_mean, 6) if rows else None,
        'read_count_ratio_cn_mean': round(read_mean, 6) if rows else None,
    },
    'host_baseline': 'autosomal_depth',
    'notes': 'Dual-method CN outputs for downstream confidence scoring in phase 2 wave 2.'
}
pathlib.Path('ebv_cn_cohort_summary.json').write_text(json.dumps(summary, indent=2))
