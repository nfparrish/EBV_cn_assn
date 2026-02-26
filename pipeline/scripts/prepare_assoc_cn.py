import csv
import json
import pathlib
import statistics
import sys

cn_table = pathlib.Path(sys.argv[1])
covariate_table = pathlib.Path(sys.argv[2])
cn_column = sys.argv[3] if len(sys.argv) > 3 else 'depth_ratio_cn_norm'

cn_rows = {}
with cn_table.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    header = reader.fieldnames or []
    if cn_column not in header:
        raise ValueError(f"CN column '{cn_column}' not found in {cn_table}")
    for row in reader:
        sid = (row.get('sample_id') or '').strip()
        if sid:
            cn_rows[sid] = row

cov_rows = {}
with covariate_table.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        sid = (row.get('sample_id') or '').strip()
        if sid:
            cov_rows[sid] = row

samples = sorted(set(cn_rows.keys()).intersection(cov_rows.keys()))
if not samples:
    raise RuntimeError('No overlapping samples between CN table and covariates')

cn_values = [float(cn_rows[sid][cn_column]) for sid in samples]
cn_mean = statistics.mean(cn_values)
cn_std = statistics.pstdev(cn_values)

def zscore(value):
    if cn_std == 0:
        return 0.0
    return (value - cn_mean) / cn_std

header = [
    'sample_id', 'cn_value', 'cn_zscore',
    'dataset', 'population', 'sex',
    'pc1', 'pc2', 'pc3', 'pc4', 'pc5',
    'sequence_lab', 'batch', 'lab_batch',
]

with pathlib.Path('association_input.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for sid in samples:
        cn_value = float(cn_rows[sid][cn_column])
        cov = cov_rows[sid]
        writer.writerow({
            'sample_id': sid,
            'cn_value': f"{cn_value:.8f}",
            'cn_zscore': f"{zscore(cn_value):.8f}",
            'dataset': cov.get('dataset', ''),
            'population': cov.get('population', ''),
            'sex': cov.get('sex', ''),
            'pc1': cov.get('pc1', ''),
            'pc2': cov.get('pc2', ''),
            'pc3': cov.get('pc3', ''),
            'pc4': cov.get('pc4', ''),
            'pc5': cov.get('pc5', ''),
            'sequence_lab': cov.get('sequence_lab', ''),
            'batch': cov.get('batch', ''),
            'lab_batch': cov.get('lab_batch', ''),
        })

report = {
    'cn_source_column': cn_column,
    'sample_count_overlap': len(samples),
    'cn_mean': cn_mean,
    'cn_std_population': cn_std,
    'cn_zscore_note': 'cn_zscore uses population std; if std is 0, all z-scores are 0',
}
pathlib.Path('association_cn_prepare_report.json').write_text(json.dumps(report, indent=2))
