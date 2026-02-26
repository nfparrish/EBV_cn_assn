import json
import pathlib
import sys

filter_counts = json.loads(pathlib.Path(sys.argv[1]).read_text())
filter_rules = json.loads(pathlib.Path(sys.argv[2]).read_text())
ebv_matrix_summary = json.loads(pathlib.Path(sys.argv[3]).read_text())
host_ingest_report = json.loads(pathlib.Path(sys.argv[4]).read_text())
host_harmonize_report = json.loads(pathlib.Path(sys.argv[5]).read_text())

report = {
    'filter_profile': filter_rules.get('profile'),
    'thresholds': filter_rules,
    'ebv': {
        'pre_filter_variant_rows': filter_counts['ebv']['unfiltered_count'],
        'post_filter_variant_rows': filter_counts['ebv']['retained_count'],
        'excluded_variant_rows': filter_counts['ebv']['excluded_count'],
        'matrix_sample_count': ebv_matrix_summary.get('sample_count'),
        'matrix_variant_count': ebv_matrix_summary.get('variant_count'),
    },
    'host': {
        'ingest_input_sample_count': host_ingest_report.get('input_sample_count'),
        'ingest_overlap_sample_count': host_ingest_report.get('overlap_sample_count'),
        'harmonized_row_count': host_harmonize_report.get('harmonized_row_count'),
        'post_filter_variant_rows': filter_counts['host']['retained_count'],
        'excluded_variant_rows': filter_counts['host']['excluded_count'],
    },
}
pathlib.Path('variant_provenance_summary.json').write_text(json.dumps(report, indent=2))

with pathlib.Path('variant_provenance_counts.tsv').open('w') as handle:
    handle.write('pathway\tstage\tcount\n')
    handle.write(f"ebv\tpre_filter\t{filter_counts['ebv']['unfiltered_count']}\n")
    handle.write(f"ebv\tpost_filter\t{filter_counts['ebv']['retained_count']}\n")
    handle.write(f"ebv\texcluded\t{filter_counts['ebv']['excluded_count']}\n")
    handle.write(f"host\tingest_input_samples\t{host_ingest_report.get('input_sample_count', 0)}\n")
    handle.write(f"host\tingest_overlap_samples\t{host_ingest_report.get('overlap_sample_count', 0)}\n")
    handle.write(f"host\tharmonized_rows\t{host_harmonize_report.get('harmonized_row_count', 0)}\n")
    handle.write(f"host\tpost_filter\t{filter_counts['host']['retained_count']}\n")
    handle.write(f"host\texcluded\t{filter_counts['host']['excluded_count']}\n")

lines = [
    '# Variant Filtering Provenance Report',
    '',
    f"Filter profile: {filter_rules.get('profile')}",
    '',
    '## EBV pathway',
    f"- Pre-filter rows: {filter_counts['ebv']['unfiltered_count']}",
    f"- Post-filter rows: {filter_counts['ebv']['retained_count']}",
    f"- Excluded rows: {filter_counts['ebv']['excluded_count']}",
    '',
    '## Host pathway',
    f"- Input samples: {host_ingest_report.get('input_sample_count')}",
    f"- Overlap samples: {host_ingest_report.get('overlap_sample_count')}",
    f"- Harmonized rows: {host_harmonize_report.get('harmonized_row_count')}",
    f"- Post-filter rows: {filter_counts['host']['retained_count']}",
    f"- Excluded rows: {filter_counts['host']['excluded_count']}",
    '',
    'This report records retained/excluded counts and thresholds for VAR-03 sign-off.',
]
pathlib.Path('variant_provenance_report.md').write_text('\n'.join(lines) + '\n')
