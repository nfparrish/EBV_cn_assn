import json
import math
import sys
from pathlib import Path

depth_path = Path(sys.argv[1])
payload = json.loads(depth_path.read_text())
metrics = payload['metrics']

host_depth = max(metrics['host_autosomal_depth'], 1e-6)
ebv_depth = max(metrics['ebv_depth'], 0.0)
total_reads = max(metrics['total_reads'], 1)
ebv_reads = max(metrics['ebv_mapped_reads'], 0)

depth_ratio_cn = round((ebv_depth / host_depth) * 2.0 * 1000.0, 6)
read_count_ratio_cn = round((ebv_reads / total_reads) * 1_000_000.0, 6)

cn_raw = {
    'sample_id': payload['sample_id'],
    'dataset': payload.get('dataset', ''),
    'population': payload.get('population', ''),
    'sex': payload.get('sex', ''),
    'sequence_lab': payload.get('sequence_lab', ''),
    'reference_accession': payload.get('reference_accession'),
    'methods': {
        'depth_ratio_cn': depth_ratio_cn,
        'read_count_ratio_cn': read_count_ratio_cn,
    },
    'formula_trace': {
        'depth_ratio': '((ebv_depth / host_autosomal_depth) * 2 * 1000)',
        'read_count_ratio': '((ebv_mapped_reads / total_reads) * 1e6)',
    },
    'support_metrics': {
        'host_autosomal_depth': host_depth,
        'ebv_depth': ebv_depth,
        'total_reads': total_reads,
        'ebv_mapped_reads': ebv_reads,
    }
}

Path(f"{payload['sample_id']}.ebv_cn_raw.json").write_text(json.dumps(cn_raw, indent=2))
