import glob
import json
import pathlib
import statistics

rows = []
for file_path in sorted(glob.glob('*.qc.json')):
    rows.append(json.loads(pathlib.Path(file_path).read_text()))

header = [
    'sample_id',
    'status',
    'dataset',
    'population',
    'sequence_lab',
    'read_count',
    'q30_rate',
    'mapping_rate',
    'mean_coverage',
    'warnings',
]

with pathlib.Path('qc_cohort_metrics.tsv').open('w') as handle:
    handle.write('\t'.join(header) + '\n')
    for row in rows:
        metric = row['metrics']
        warning_str = ','.join(row.get('warnings', []))
        values = [
            row.get('sample_id', ''),
            row.get('status', ''),
            row.get('dataset', ''),
            row.get('population', ''),
            row.get('sequence_lab', ''),
            str(metric.get('read_count', '')),
            str(metric.get('q30_rate', '')),
            str(metric.get('mapping_rate', '')),
            str(metric.get('mean_coverage', '')),
            warning_str,
        ]
        handle.write('\t'.join(values) + '\n')

q30_values = [r['metrics']['q30_rate'] for r in rows] if rows else []
map_values = [r['metrics']['mapping_rate'] for r in rows] if rows else []
cov_values = [r['metrics']['mean_coverage'] for r in rows] if rows else []

summary = {
    'sample_count': len(rows),
    'status_counts': {
        'processed': sum(1 for r in rows if r.get('status') == 'processed'),
        'skipped': sum(1 for r in rows if r.get('status') == 'skipped'),
    },
    'metrics': {
        'q30_rate_mean': round(statistics.mean(q30_values), 4) if q30_values else None,
        'mapping_rate_mean': round(statistics.mean(map_values), 4) if map_values else None,
        'mean_coverage_mean': round(statistics.mean(cov_values), 3) if cov_values else None,
    },
    'warning_counts': {
        'low_coverage': sum(1 for r in rows if 'low_coverage' in r.get('warnings', [])),
        'low_q30': sum(1 for r in rows if 'low_q30' in r.get('warnings', [])),
    },
}

pathlib.Path('qc_cohort_summary.json').write_text(json.dumps(summary, indent=2))

report_lines = [
    '# Cohort QC Report',
    '',
    f"Samples analyzed: {summary['sample_count']}",
    f"Processed: {summary['status_counts']['processed']}",
    f"Skipped: {summary['status_counts']['skipped']}",
    '',
    '## Mean Metrics',
    f"- Mean Q30 rate: {summary['metrics']['q30_rate_mean']}",
    f"- Mean mapping rate: {summary['metrics']['mapping_rate_mean']}",
    f"- Mean coverage: {summary['metrics']['mean_coverage_mean']}",
    '',
    '## Warning Counts',
    f"- low_coverage: {summary['warning_counts']['low_coverage']}",
    f"- low_q30: {summary['warning_counts']['low_q30']}",
    '',
    'This report is intended for pre-CN intake QC validation.',
]
pathlib.Path('qc_cohort_report.md').write_text('\n'.join(report_lines) + '\n')
