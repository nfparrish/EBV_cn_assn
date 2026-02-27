import csv
import json
import pathlib
import sys


def resolve_path(path_value, project_dir):
    if not path_value:
        return None
    candidate = pathlib.Path(path_value)
    if candidate.is_absolute():
        return candidate
    if candidate.exists():
        return candidate
    return pathlib.Path(project_dir) / candidate


manifest_arg = sys.argv[1]
covariates_arg = sys.argv[2] if len(sys.argv) > 2 else ''
project_dir = sys.argv[3] if len(sys.argv) > 3 else '.'

manifest_path = resolve_path(manifest_arg, project_dir)
if manifest_path is None or not manifest_path.exists():
    raise FileNotFoundError(f"Manifest path not found: {manifest_arg}")

covariates_path = resolve_path(covariates_arg, project_dir) if covariates_arg else None
if covariates_path and not covariates_path.exists():
    raise FileNotFoundError(f"Covariates path not found: {covariates_arg}")

manifest_rows = []
with manifest_path.open() as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        sid = (row.get('sample_id') or '').strip()
        if sid:
            manifest_rows.append(row)

manifest_by_sample = {r['sample_id'].strip(): r for r in manifest_rows}
manifest_samples = sorted(manifest_by_sample.keys())

source_rows = {}
source_mode = 'input_table'

if not covariates_path:
    raise ValueError(
        'Empirical association covariates are required. '
        'Provide --association_covariates_input with sample-level covariates including pc1..pc5.'
    )

with covariates_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        sid = (row.get('sample_id') or '').strip()
        if sid:
            source_rows[sid] = row

missing_covariate_samples = [sid for sid in manifest_samples if sid not in source_rows]
if missing_covariate_samples:
    raise ValueError(
        'Covariate table is missing manifest samples: '
        + ', '.join(missing_covariate_samples[:20])
        + (' ...' if len(missing_covariate_samples) > 20 else '')
    )

header = [
    'sample_id', 'dataset', 'population', 'sex',
    'pc1', 'pc2', 'pc3', 'pc4', 'pc5',
    'sequence_lab', 'batch', 'lab_batch',
    'covariate_mode',
]

output_rows = []
for sid in manifest_samples:
    manifest_meta = manifest_by_sample[sid]
    base = source_rows.get(sid, {})

    sex = (base.get('sex') or manifest_meta.get('sex') or 'U').strip() or 'U'
    sequence_lab = (base.get('sequence_lab') or manifest_meta.get('sequence_lab') or '').strip()
    batch = (base.get('batch') or '').strip()

    lab_batch = (base.get('lab_batch') or '').strip()
    if not lab_batch:
        left = sequence_lab if sequence_lab else ''
        right = batch if batch else ''
        lab_batch = f"{left}__{right}" if (left or right) else ''

    pc_values = {}
    for field in ['pc1', 'pc2', 'pc3', 'pc4', 'pc5']:
        raw = (base.get(field) or '').strip()
        if raw == '':
            raise ValueError(f"Missing required covariate '{field}' for sample {sid}")
        try:
            pc_values[field] = f"{float(raw):.6f}"
        except ValueError as exc:
            raise ValueError(f"Invalid numeric covariate '{field}' for sample {sid}: {raw}") from exc

    row = {
        'sample_id': sid,
        'dataset': manifest_meta.get('dataset', ''),
        'population': manifest_meta.get('population', ''),
        'sex': sex,
        'pc1': pc_values['pc1'],
        'pc2': pc_values['pc2'],
        'pc3': pc_values['pc3'],
        'pc4': pc_values['pc4'],
        'pc5': pc_values['pc5'],
        'sequence_lab': sequence_lab,
        'batch': batch,
        'lab_batch': lab_batch,
        'covariate_mode': 'empirical_input',
    }
    output_rows.append(row)

with pathlib.Path('association_covariates.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in output_rows:
        writer.writerow(row)

summary = {
    'source_mode': source_mode,
    'manifest_sample_count': len(manifest_samples),
    'covariate_sample_count': len(output_rows),
    'missing_from_covariates_count': 0,
    'fields': header,
}
pathlib.Path('association_covariates_report.json').write_text(json.dumps(summary, indent=2))
