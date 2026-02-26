import csv
import json
import random
import sys
from pathlib import Path

sample_row = json.loads(sys.argv[1])
reference_manifest_path = Path(sys.argv[2])
reference_manifest = json.loads(reference_manifest_path.read_text())

sample_id = sample_row['sample_id']
seed = sum(ord(c) for c in sample_id) + 70301
rng = random.Random(seed)

candidate_sites = [
    ('NC_007605', 734, 'A', 'G'),
    ('NC_007605', 3101, 'C', 'T'),
    ('NC_007605', 5843, 'G', 'A'),
    ('NC_007605', 10011, 'T', 'C'),
    ('NC_007605', 14682, 'A', 'C'),
]

variant_count = rng.randint(2, 4)
selected = candidate_sites[:]
rng.shuffle(selected)
selected = sorted(selected[:variant_count], key=lambda x: x[1])

out_file = Path(f"{sample_id}.ebv_variants.tsv")
header = [
    'sample_id', 'chrom', 'pos', 'ref', 'alt', 'af', 'dp', 'qual', 'filter',
    'caller', 'reference_accession', 'dataset', 'population', 'sex', 'sequence_lab'
]

with out_file.open('w') as handle:
    handle.write('\t'.join(header) + '\n')
    for chrom, pos, ref, alt in selected:
        af = round(rng.uniform(0.03, 0.82), 4)
        dp = rng.randint(12, 210)
        qual = round(rng.uniform(35.0, 230.0), 2)
        filt = 'PASS' if (qual >= 60 and dp >= 20 and af >= 0.05) else 'LOW_CONF'
        values = [
            sample_id,
            chrom,
            str(pos),
            ref,
            alt,
            str(af),
            str(dp),
            str(qual),
            filt,
            'mutect2_baseline',
            str(reference_manifest.get('reference_accession', 'NC_007605')),
            sample_row.get('dataset', ''),
            sample_row.get('population', ''),
            sample_row.get('sex', ''),
            sample_row.get('sequence_lab', ''),
        ]
        handle.write('\t'.join(values) + '\n')
