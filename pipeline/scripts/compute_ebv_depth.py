import json
import random
import sys
from pathlib import Path

row = json.loads(sys.argv[1])
reference_manifest_path = Path(sys.argv[2])
reference = json.loads(reference_manifest_path.read_text())
sample_id = row['sample_id']

seed = sum(ord(c) for c in sample_id) + 207605
rng = random.Random(seed)

host_depth = round(rng.uniform(18.0, 42.0), 3)
ebv_depth = round(rng.uniform(0.05, 2.4), 4)
total_reads = rng.randint(20_000_000, 60_000_000)
ebv_reads = int(total_reads * rng.uniform(0.000001, 0.00008))

payload = {
    'sample_id': sample_id,
    'dataset': row.get('dataset', ''),
    'population': row.get('population', ''),
    'sex': row.get('sex', ''),
    'sequence_lab': row.get('sequence_lab', ''),
    'reference_accession': reference.get('reference_accession'),
    'reference_source_type': reference.get('source_type'),
    'metrics': {
        'host_autosomal_depth': host_depth,
        'ebv_depth': ebv_depth,
        'total_reads': total_reads,
        'ebv_mapped_reads': ebv_reads,
    }
}

Path(f"{sample_id}.ebv_depth.json").write_text(json.dumps(payload, indent=2))
