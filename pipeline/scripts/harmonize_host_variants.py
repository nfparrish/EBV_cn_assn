import csv
import json
import pathlib
import sys

manifest_path = pathlib.Path(sys.argv[1])
staged_variants_path = pathlib.Path(sys.argv[2])

manifest_rows = []
with manifest_path.open() as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        sample_id = row.get('sample_id', '').strip()
        if sample_id:
            manifest_rows.append(row)

manifest_by_sample = {row['sample_id'].strip(): row for row in manifest_rows}
manifest_sample_set = set(manifest_by_sample.keys())

staged_rows = []
with staged_variants_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        staged_rows.append(row)

input_sample_set = {row.get('sample_id', '').strip() for row in staged_rows if row.get('sample_id', '').strip()}
overlap_samples = sorted(input_sample_set.intersection(manifest_sample_set))
missing_manifest_samples = sorted(input_sample_set - manifest_sample_set)
missing_input_samples = sorted(manifest_sample_set - input_sample_set)

out_header = [
    'sample_id', 'dataset', 'population', 'sex', 'sequence_lab',
    'variant_id', 'chrom', 'pos', 'ref', 'alt', 'gt', 'dosage'
]

harmonized_rows = 0
with pathlib.Path('host_variant_harmonized.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=out_header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in staged_rows:
        sample_id = row.get('sample_id', '').strip()
        if sample_id not in manifest_by_sample:
            continue
        meta = manifest_by_sample[sample_id]
        writer.writerow({
            'sample_id': sample_id,
            'dataset': meta.get('dataset', ''),
            'population': meta.get('population', ''),
            'sex': meta.get('sex', ''),
            'sequence_lab': meta.get('sequence_lab', ''),
            'variant_id': row.get('variant_id', ''),
            'chrom': row.get('chrom', ''),
            'pos': row.get('pos', ''),
            'ref': row.get('ref', ''),
            'alt': row.get('alt', ''),
            'gt': row.get('gt', ''),
            'dosage': row.get('dosage', ''),
        })
        harmonized_rows += 1

with pathlib.Path('host_variant_sample_mismatch.tsv').open('w') as handle:
    handle.write('category\tsample_id\n')
    for sample_id in missing_manifest_samples:
        handle.write(f"input_not_in_manifest\t{sample_id}\n")
    for sample_id in missing_input_samples:
        handle.write(f"manifest_not_in_input\t{sample_id}\n")

summary = {
    'manifest_sample_count': len(manifest_sample_set),
    'input_sample_count': len(input_sample_set),
    'overlap_sample_count': len(overlap_samples),
    'harmonized_row_count': harmonized_rows,
    'missing_from_manifest_count': len(missing_manifest_samples),
    'missing_from_input_count': len(missing_input_samples),
}
pathlib.Path('host_variant_harmonize_report.json').write_text(json.dumps(summary, indent=2))
