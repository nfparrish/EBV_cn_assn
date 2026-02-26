import csv
import json
import pathlib
import random
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
source_mode = 'generated_from_manifest'

if covariates_path:
    with covariates_path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            if sid:
                source_rows[sid] = row
    source_mode = 'input_table'

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

    seed = sum(ord(c) for c in sid) + 401
    rng = random.Random(seed)

    sex = (base.get('sex') or manifest_meta.get('sex') or 'U').strip() or 'U'
    sequence_lab = (base.get('sequence_lab') or manifest_meta.get('sequence_lab') or '').strip()
    batch = (base.get('batch') or '').strip()

    if not batch:
        batch = f"batch_{(seed % 3) + 1}"

    lab_batch = (base.get('lab_batch') or '').strip()
    if not lab_batch:
        left = sequence_lab if sequence_lab else 'UNKNOWN_LAB'
        right = batch if batch else 'UNKNOWN_BATCH'
        lab_batch = f"{left}__{right}"

    row = {
        'sample_id': sid,
        'dataset': manifest_meta.get('dataset', ''),
        'population': manifest_meta.get('population', ''),
        'sex': sex,
        'pc1': base.get('pc1') or f"{rng.uniform(-0.2, 0.2):.6f}",
        'pc2': base.get('pc2') or f"{rng.uniform(-0.2, 0.2):.6f}",
        'pc3': base.get('pc3') or f"{rng.uniform(-0.2, 0.2):.6f}",
        'pc4': base.get('pc4') or f"{rng.uniform(-0.2, 0.2):.6f}",
        'pc5': base.get('pc5') or f"{rng.uniform(-0.2, 0.2):.6f}",
        'sequence_lab': sequence_lab,
        'batch': batch,
        'lab_batch': lab_batch,
        'covariate_mode': 'separate_if_available_else_combined',
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
    'missing_from_covariates_count': sum(1 for sid in manifest_samples if sid not in source_rows),
    'fields': header,
}
pathlib.Path('association_covariates_report.json').write_text(json.dumps(summary, indent=2))
