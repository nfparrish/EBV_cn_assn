import glob
import json
import pathlib

rows = []
for file_path in glob.glob('*.status.json'):
    rows.append(json.loads(pathlib.Path(file_path).read_text()))

rows = sorted(rows, key=lambda item: item['sample_id'])

header = ['sample_id', 'status', 'reason', 'dataset', 'population', 'sequence_lab']
with pathlib.Path('sample_status_registry.tsv').open('w') as handle:
    handle.write('\t'.join(header) + '\n')
    for row in rows:
        handle.write('\t'.join(str(row.get(col, '')) for col in header) + '\n')

summary = {
    'total': len(rows),
    'processed': sum(1 for row in rows if row.get('status') == 'processed'),
    'skipped': sum(1 for row in rows if row.get('status') == 'skipped'),
}
pathlib.Path('sample_status_summary.json').write_text(json.dumps(summary, indent=2))
