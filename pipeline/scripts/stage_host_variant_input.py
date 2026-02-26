import csv
import json
import pathlib
import sys

manifest_path = pathlib.Path(sys.argv[1])
host_input_path = pathlib.Path(sys.argv[2])
host_variant_format = sys.argv[3].lower()

if host_variant_format not in {'tsv', 'vcf', 'plink'}:
    raise ValueError(f"Unsupported host_variant_format: {host_variant_format}")

manifest_samples = []
with manifest_path.open() as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        sample_id = row.get('sample_id', '').strip()
        if sample_id:
            manifest_samples.append(sample_id)

manifest_sample_set = set(manifest_samples)

rows = []
input_header = []
with host_input_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    input_header = reader.fieldnames or []
    required = {'sample_id', 'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gt', 'dosage'}
    missing = [col for col in sorted(required) if col not in input_header]
    if missing:
        raise ValueError(f"Host input missing required columns: {', '.join(missing)}")
    for row in reader:
        rows.append(row)

with pathlib.Path('host_variant_staged.tsv').open('w') as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=['sample_id', 'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gt', 'dosage'],
        delimiter='\t',
        lineterminator='\n',
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({
            'sample_id': row.get('sample_id', ''),
            'variant_id': row.get('variant_id', ''),
            'chrom': row.get('chrom', ''),
            'pos': row.get('pos', ''),
            'ref': row.get('ref', ''),
            'alt': row.get('alt', ''),
            'gt': row.get('gt', ''),
            'dosage': row.get('dosage', ''),
        })

input_sample_set = {row.get('sample_id', '').strip() for row in rows if row.get('sample_id', '').strip()}
in_manifest = sorted(input_sample_set.intersection(manifest_sample_set))
missing_from_manifest = sorted(input_sample_set - manifest_sample_set)
missing_from_input = sorted(manifest_sample_set - input_sample_set)

summary = {
    'host_variant_format': host_variant_format,
    'input_path': str(host_input_path),
    'input_columns': input_header,
    'manifest_sample_count': len(manifest_sample_set),
    'input_sample_count': len(input_sample_set),
    'overlap_sample_count': len(in_manifest),
    'missing_from_manifest_count': len(missing_from_manifest),
    'missing_from_input_count': len(missing_from_input),
}
pathlib.Path('host_variant_ingest_report.json').write_text(json.dumps(summary, indent=2))

with pathlib.Path('host_variant_sample_map.tsv').open('w') as handle:
    handle.write('sample_id\tin_manifest\n')
    for sample_id in sorted(input_sample_set):
        handle.write(f"{sample_id}\t{1 if sample_id in manifest_sample_set else 0}\n")
