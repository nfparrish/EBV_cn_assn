import csv
import json
import pathlib
import statistics
import sys

rows = []
confidence_input = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path('ebv_cn_confidence.tsv')
with confidence_input.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        rows.append(row)

scores = [int(row.get('confidence_score', 0)) for row in rows] if rows else []
low_rows = [row for row in rows if row.get('confidence_tier') == 'low']

with pathlib.Path('ebv_cn_low_confidence.tsv').open('w') as handle:
    header = [
        'sample_id', 'dataset', 'population', 'sex', 'sequence_lab',
        'confidence_score', 'confidence_tier', 'confidence_flags'
    ]
    handle.write('\t'.join(header) + '\n')
    for row in low_rows:
        values = [
            row.get('sample_id', ''),
            row.get('dataset', ''),
            row.get('population', ''),
            row.get('sex', ''),
            row.get('sequence_lab', ''),
            row.get('confidence_score', ''),
            row.get('confidence_tier', ''),
            row.get('confidence_flags', ''),
        ]
        handle.write('\t'.join(values) + '\n')

qc_summary = {
    'sample_count': len(rows),
    'confidence_score_mean': round(statistics.mean(scores), 3) if scores else None,
    'confidence_score_median': round(statistics.median(scores), 3) if scores else None,
    'low_confidence_count': len(low_rows),
    'low_confidence_rate': round(len(low_rows) / len(rows), 4) if rows else None,
    'tier_counts': {
        'high': sum(1 for r in rows if r.get('confidence_tier') == 'high'),
        'medium': sum(1 for r in rows if r.get('confidence_tier') == 'medium'),
        'low': sum(1 for r in rows if r.get('confidence_tier') == 'low'),
    },
    'flag_counts': {
        'low_depth_signal': sum(1 for r in rows if 'low_depth_signal' in r.get('confidence_flags', '')),
        'low_readcount_signal': sum(1 for r in rows if 'low_readcount_signal' in r.get('confidence_flags', '')),
        'method_discordance_moderate': sum(1 for r in rows if 'method_discordance_moderate' in r.get('confidence_flags', '')),
        'method_discordance_high': sum(1 for r in rows if 'method_discordance_high' in r.get('confidence_flags', '')),
        'extreme_normalized_signal': sum(1 for r in rows if 'extreme_normalized_signal' in r.get('confidence_flags', '')),
    },
}

pathlib.Path('ebv_cn_qc_summary.json').write_text(json.dumps(qc_summary, indent=2))

report_lines = [
    '# EBV CN QC Report',
    '',
    f"Samples analyzed: {qc_summary['sample_count']}",
    f"Low-confidence calls: {qc_summary['low_confidence_count']}",
    f"Low-confidence rate: {qc_summary['low_confidence_rate']}",
    '',
    '## Confidence Score',
    f"- Mean: {qc_summary['confidence_score_mean']}",
    f"- Median: {qc_summary['confidence_score_median']}",
    '',
    '## Tier Counts',
    f"- high: {qc_summary['tier_counts']['high']}",
    f"- medium: {qc_summary['tier_counts']['medium']}",
    f"- low: {qc_summary['tier_counts']['low']}",
    '',
    '## Flag Counts',
    f"- low_depth_signal: {qc_summary['flag_counts']['low_depth_signal']}",
    f"- low_readcount_signal: {qc_summary['flag_counts']['low_readcount_signal']}",
    f"- method_discordance_moderate: {qc_summary['flag_counts']['method_discordance_moderate']}",
    f"- method_discordance_high: {qc_summary['flag_counts']['method_discordance_high']}",
    f"- extreme_normalized_signal: {qc_summary['flag_counts']['extreme_normalized_signal']}",
    '',
    'Pilot QC summary is intended for phase sign-off and threshold tuning.',
]
pathlib.Path('ebv_cn_qc_report.md').write_text('\n'.join(report_lines) + '\n')
