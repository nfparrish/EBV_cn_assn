import csv
import glob
import json
import pathlib
import sys


def as_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


ebv_min_af = as_float(sys.argv[1], 0.05)
ebv_min_dp = as_float(sys.argv[2], 20)
ebv_min_qual = as_float(sys.argv[3], 60)
host_min_dosage = as_float(sys.argv[4], 1.0)
filter_profile = sys.argv[5] if len(sys.argv) > 5 else 'balanced'
host_harmonized_path = pathlib.Path(sys.argv[6])

rules = {
    'profile': filter_profile,
    'ebv': {
        'min_af': ebv_min_af,
        'min_dp': ebv_min_dp,
        'min_qual': ebv_min_qual,
    },
    'host': {
        'min_dosage': host_min_dosage,
    },
}
pathlib.Path('variant_filter_rules.json').write_text(json.dumps(rules, indent=2))

ebv_files = sorted(glob.glob('*.ebv_variants.tsv'))
ebv_rows = []
for file_path in ebv_files:
    with pathlib.Path(file_path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            ebv_rows.append(row)

if not ebv_rows:
    raise RuntimeError('No EBV variant rows found for filtering')

ebv_header = [
    'sample_id', 'chrom', 'pos', 'ref', 'alt', 'af', 'dp', 'qual', 'filter',
    'caller', 'reference_accession', 'dataset', 'population', 'sex', 'sequence_lab'
]
with pathlib.Path('ebv_variants_unfiltered.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=ebv_header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in ebv_rows:
        writer.writerow({k: row.get(k, '') for k in ebv_header})

ebv_retained = []
ebv_excluded = []
for row in ebv_rows:
    af = as_float(row.get('af'), 0.0)
    dp = as_float(row.get('dp'), 0.0)
    qual = as_float(row.get('qual'), 0.0)
    if af >= ebv_min_af and dp >= ebv_min_dp and qual >= ebv_min_qual:
        ebv_retained.append(row)
    else:
        ebv_excluded.append(row)

with pathlib.Path('ebv_variants_filtered.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=ebv_header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in ebv_retained:
        writer.writerow({k: row.get(k, '') for k in ebv_header})

host_rows = []
with host_harmonized_path.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    host_header = reader.fieldnames or []
    for row in reader:
        host_rows.append(row)

if not host_header:
    raise RuntimeError('Host harmonized input missing header')

host_retained = []
host_excluded = []
for row in host_rows:
    dosage = as_float(row.get('dosage'), 0.0)
    if dosage >= host_min_dosage:
        host_retained.append(row)
    else:
        host_excluded.append(row)

with pathlib.Path('host_variant_filtered.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=host_header, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in host_retained:
        writer.writerow({k: row.get(k, '') for k in host_header})

counts = {
    'ebv': {
        'unfiltered_count': len(ebv_rows),
        'retained_count': len(ebv_retained),
        'excluded_count': len(ebv_excluded),
        'rules': {
            'min_af': ebv_min_af,
            'min_dp': ebv_min_dp,
            'min_qual': ebv_min_qual,
        },
    },
    'host': {
        'unfiltered_count': len(host_rows),
        'retained_count': len(host_retained),
        'excluded_count': len(host_excluded),
        'rules': {
            'min_dosage': host_min_dosage,
        },
    },
}
pathlib.Path('variant_filter_counts.json').write_text(json.dumps(counts, indent=2))
