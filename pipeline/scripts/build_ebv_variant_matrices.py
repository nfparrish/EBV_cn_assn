import csv
import glob
import json
import pathlib
import statistics
import sys

files = sorted(glob.glob('*.ebv_variants.tsv'))
rows = []
samples = set()
variants = set()

for file_path in files:
    with pathlib.Path(file_path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            rows.append(row)
            samples.add(row['sample_id'])
            variant_id = f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
            variants.add(variant_id)

samples = sorted(samples)
variants = sorted(variants, key=lambda v: (v.split(':')[0], int(v.split(':')[1]), v))

af_threshold = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05

value_map = {}
for row in rows:
    variant_id = f"{row['chrom']}:{row['pos']}:{row['ref']}>{row['alt']}"
    key = (row['sample_id'], variant_id)
    af = float(row['af'])
    if key not in value_map or af > value_map[key]['af']:
        value_map[key] = {
            'af': af,
            'dp': int(row['dp']),
            'qual': float(row['qual']),
            'filter': row['filter'],
        }

with pathlib.Path('ebv_variant_catalog.tsv').open('w') as handle:
    catalog_header = ['variant_id', 'chrom', 'pos', 'ref', 'alt', 'carrier_count', 'mean_af']
    handle.write('\t'.join(catalog_header) + '\n')
    for variant_id in variants:
        chrom, pos, change = variant_id.split(':', 2)
        ref, alt = change.split('>')
        carrier_afs = [
            value_map[(sample_id, variant_id)]['af']
            for sample_id in samples
            if (sample_id, variant_id) in value_map
        ]
        carrier_count = len([x for x in carrier_afs if x >= af_threshold])
        mean_af = round(statistics.mean(carrier_afs), 6) if carrier_afs else 0.0
        handle.write('\t'.join([
            variant_id,
            chrom,
            pos,
            ref,
            alt,
            str(carrier_count),
            str(mean_af),
        ]) + '\n')

af_header = ['sample_id'] + variants
with pathlib.Path('ebv_variant_af_matrix.tsv').open('w') as handle:
    handle.write('\t'.join(af_header) + '\n')
    for sample_id in samples:
        values = [sample_id]
        for variant_id in variants:
            af = value_map.get((sample_id, variant_id), {}).get('af', 0.0)
            values.append(str(round(af, 6)))
        handle.write('\t'.join(values) + '\n')

binary_header = ['sample_id'] + variants
with pathlib.Path('ebv_variant_binary_matrix.tsv').open('w') as handle:
    handle.write('\t'.join(binary_header) + '\n')
    for sample_id in samples:
        values = [sample_id]
        for variant_id in variants:
            af = value_map.get((sample_id, variant_id), {}).get('af', 0.0)
            values.append('1' if af >= af_threshold else '0')
        handle.write('\t'.join(values) + '\n')

summary = {
    'caller': 'mutect2_baseline',
    'sample_count': len(samples),
    'variant_count': len(variants),
    'binary_af_threshold': af_threshold,
    'filter_counts': {
        'PASS': sum(1 for row in rows if row.get('filter') == 'PASS'),
        'LOW_CONF': sum(1 for row in rows if row.get('filter') == 'LOW_CONF'),
    },
    'notes': 'Phase 3 Wave 1 dual-matrix EBV variant aggregation output.',
}
pathlib.Path('ebv_variant_matrix_summary.json').write_text(json.dumps(summary, indent=2))
