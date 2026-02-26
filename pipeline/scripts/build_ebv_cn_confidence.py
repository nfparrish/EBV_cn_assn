import csv
import json
import pathlib
import sys


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def score_row(row):
    depth_raw = to_float(row.get('depth_ratio_cn'))
    read_raw = to_float(row.get('read_count_ratio_cn'))
    depth_norm = to_float(row.get('depth_ratio_cn_norm'))
    read_norm = to_float(row.get('read_count_ratio_cn_norm'))

    flags = []
    score = 100

    if depth_raw < 0.5:
        score -= 20
        flags.append('low_depth_signal')

    if read_raw < 1.0:
        score -= 20
        flags.append('low_readcount_signal')

    discordance = abs(depth_norm - read_norm)
    if discordance > 0.5:
        score -= 30
        flags.append('method_discordance_high')
    elif discordance > 0.25:
        score -= 15
        flags.append('method_discordance_moderate')

    if depth_norm < 0.25 or read_norm < 0.25:
        score -= 20
        flags.append('extreme_normalized_signal')

    score = max(0, min(100, score))
    if score >= 80:
        tier = 'high'
    elif score >= 55:
        tier = 'medium'
    else:
        tier = 'low'

    return score, tier, flags


normalized_input = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path('ebv_cn_normalized.tsv')
summary_input = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else pathlib.Path('ebv_cn_cohort_summary.json')

normalized_rows = []
with normalized_input.open() as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        normalized_rows.append(row)

cohort_summary = json.loads(summary_input.read_text())

header = [
    'sample_id', 'dataset', 'population', 'sex', 'sequence_lab', 'reference_accession',
    'depth_ratio_cn', 'read_count_ratio_cn', 'depth_ratio_cn_norm', 'read_count_ratio_cn_norm',
    'confidence_score', 'confidence_tier', 'confidence_flags'
]

low_confidence_count = 0
medium_confidence_count = 0
high_confidence_count = 0

with pathlib.Path('ebv_cn_confidence.tsv').open('w') as handle:
    handle.write('\t'.join(header) + '\n')
    for row in normalized_rows:
        score, tier, flags = score_row(row)

        if tier == 'high':
            high_confidence_count += 1
        elif tier == 'medium':
            medium_confidence_count += 1
        else:
            low_confidence_count += 1

        values = [
            row.get('sample_id', ''),
            row.get('dataset', ''),
            row.get('population', ''),
            row.get('sex', ''),
            row.get('sequence_lab', ''),
            row.get('reference_accession', ''),
            row.get('depth_ratio_cn', ''),
            row.get('read_count_ratio_cn', ''),
            row.get('depth_ratio_cn_norm', ''),
            row.get('read_count_ratio_cn_norm', ''),
            str(score),
            tier,
            ','.join(flags),
        ]
        handle.write('\t'.join(values) + '\n')

confidence_summary = {
    'sample_count': len(normalized_rows),
    'confidence_counts': {
        'high': high_confidence_count,
        'medium': medium_confidence_count,
        'low': low_confidence_count,
    },
    'rules': {
        'score_start': 100,
        'depth_raw_lt_0_5_penalty': 20,
        'read_raw_lt_1_penalty': 20,
        'discordance_gt_0_5_penalty': 30,
        'discordance_gt_0_25_penalty': 15,
        'extreme_norm_lt_0_25_penalty': 20,
        'tier_high_min': 80,
        'tier_medium_min': 55,
    },
    'cohort_means': cohort_summary.get('methods', {}),
}

pathlib.Path('ebv_cn_confidence_summary.json').write_text(json.dumps(confidence_summary, indent=2))
