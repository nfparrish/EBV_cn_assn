import csv
import glob
import json
import pathlib
import sys


def to_float(value, default=0.0):
    try:
        return float(value)
    except Exception:
        return default


def resolve_path(path_value, project_dir):
    if not path_value:
        return None
    candidate = pathlib.Path(path_value)
    if candidate.is_absolute() and candidate.exists():
        return candidate
    if candidate.exists():
        return candidate
    rel = pathlib.Path(project_dir) / candidate
    return rel if rel.exists() else None


def parse_expected_coding_variants(raw):
    items = []
    for token in (raw or '').split(','):
        value = token.strip()
        if value:
            items.append(value)
    return sorted(set(items))


def load_sample_map(path):
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return sorted([(row.get('sample_id') or '').strip() for row in reader if (row.get('sample_id') or '').strip()])


def load_observed_variants(min_af, min_dp, min_qual):
    observed = {}

    for file_path in sorted(glob.glob('*.ebv_variants.tsv')):
        with pathlib.Path(file_path).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                sample_id = (row.get('sample_id') or '').strip()
                if not sample_id:
                    continue

                chrom = (row.get('chrom') or '').strip()
                pos = (row.get('pos') or '').strip()
                ref = (row.get('ref') or '').strip()
                alt = (row.get('alt') or '').strip()
                if not chrom or not pos or not ref or not alt:
                    continue

                variant_id = f"{chrom}:{pos}:{ref}>{alt}"
                af = to_float(row.get('af'), 0.0)
                dp = to_float(row.get('dp'), 0.0)
                qual = to_float(row.get('qual'), 0.0)
                status = (row.get('filter') or '').strip().upper()

                if af < min_af or dp < min_dp or qual < min_qual:
                    continue
                if status not in ('PASS', ''):
                    continue

                sample_variants = observed.setdefault(sample_id, {})
                if variant_id not in sample_variants or af > sample_variants[variant_id]['af']:
                    sample_variants[variant_id] = {'af': af, 'dp': dp, 'qual': qual}

    return observed


def load_donor_expectations(path):
    donor_expect = {}
    if not path:
        return donor_expect

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            sample_id = (row.get('sample_id') or '').strip()
            variant_id = (row.get('variant_id') or '').strip()
            if not sample_id or not variant_id:
                continue
            donor_expect.setdefault(sample_id, set()).add(variant_id)

    return donor_expect


def infer_donor_specific_expectations(observed):
    carriers = {}
    for sample_id, sample_vars in observed.items():
        for variant_id in sample_vars.keys():
            carriers.setdefault(variant_id, set()).add(sample_id)

    donor_specific = {}
    for variant_id, sample_ids in carriers.items():
        if len(sample_ids) == 1:
            sid = next(iter(sample_ids))
            donor_specific.setdefault(sid, set()).add(variant_id)

    return donor_specific


def main():
    sample_map_path = pathlib.Path(sys.argv[1])
    expected_coding_raw = sys.argv[2] if len(sys.argv) > 2 else ''
    donor_expect_arg = sys.argv[3] if len(sys.argv) > 3 else ''
    project_dir = sys.argv[4] if len(sys.argv) > 4 else '.'
    min_af = to_float(sys.argv[5], 0.05) if len(sys.argv) > 5 else 0.05
    min_dp = to_float(sys.argv[6], 20.0) if len(sys.argv) > 6 else 20.0
    min_qual = to_float(sys.argv[7], 60.0) if len(sys.argv) > 7 else 60.0

    sample_ids = load_sample_map(sample_map_path)
    expected_coding = parse_expected_coding_variants(expected_coding_raw)
    observed = load_observed_variants(min_af=min_af, min_dp=min_dp, min_qual=min_qual)

    donor_expect_file = resolve_path(donor_expect_arg, project_dir)
    donor_expect = load_donor_expectations(donor_expect_file) if donor_expect_file else {}
    donor_expect_source = 'input_table' if donor_expect_file else 'inferred_unique_variants'
    if not donor_expect:
        donor_expect = infer_donor_specific_expectations(observed)

    rows = []

    for sid in sample_ids:
        sample_variants = observed.get(sid, {})

        for variant_id in expected_coding:
            evidence = sample_variants.get(variant_id)
            rows.append({
                'sample_id': sid,
                'check_class': 'ebv_expected_coding',
                'variant_id': variant_id,
                'present': 'true' if evidence else 'false',
                'evidence_af': str(evidence['af']) if evidence else '',
                'evidence_dp': str(evidence['dp']) if evidence else '',
                'evidence_qual': str(evidence['qual']) if evidence else '',
            })

        for variant_id in sorted(donor_expect.get(sid, set())):
            evidence = sample_variants.get(variant_id)
            rows.append({
                'sample_id': sid,
                'check_class': 'donor_specific_coding',
                'variant_id': variant_id,
                'present': 'true' if evidence else 'false',
                'evidence_af': str(evidence['af']) if evidence else '',
                'evidence_dp': str(evidence['dp']) if evidence else '',
                'evidence_qual': str(evidence['qual']) if evidence else '',
            })

    with pathlib.Path('rna_variant_presence_checks.tsv').open('w') as handle:
        header = ['sample_id', 'check_class', 'variant_id', 'present', 'evidence_af', 'evidence_dp', 'evidence_qual']
        writer = csv.DictWriter(handle, fieldnames=header, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    by_sample = {}
    for row in rows:
        sid = row['sample_id']
        slot = by_sample.setdefault(sid, {
            'ebv_expected_total': 0,
            'ebv_expected_present': 0,
            'donor_expected_total': 0,
            'donor_expected_present': 0,
        })
        present = row['present'] == 'true'
        if row['check_class'] == 'ebv_expected_coding':
            slot['ebv_expected_total'] += 1
            slot['ebv_expected_present'] += 1 if present else 0
        else:
            slot['donor_expected_total'] += 1
            slot['donor_expected_present'] += 1 if present else 0

    summary = {
        'sample_count': len(sample_ids),
        'filters': {
            'min_af': min_af,
            'min_dp': min_dp,
            'min_qual': min_qual,
        },
        'expected_ebv_coding_variant_count': len(expected_coding),
        'donor_expectation_source': donor_expect_source,
        'samples': by_sample,
        'overall_flags': {
            'all_expected_ebv_coding_present': all(
                (v['ebv_expected_present'] == v['ebv_expected_total'])
                for v in by_sample.values()
            ) if by_sample else True,
            'all_expected_donor_specific_present': all(
                (v['donor_expected_present'] == v['donor_expected_total'])
                for v in by_sample.values()
            ) if by_sample else True,
        },
    }
    pathlib.Path('rna_variant_presence_summary.json').write_text(json.dumps(summary, indent=2))

    report_lines = [
        '# RNA Variant Presence QC',
        '',
        '## Scope',
        f"- Samples evaluated: {len(sample_ids)}",
        f"- Expected EBV coding variants per sample: {len(expected_coding)}",
        f"- Donor-specific expectation source: {donor_expect_source}",
        '',
        '## Global flags',
        f"- all_expected_ebv_coding_present: {summary['overall_flags']['all_expected_ebv_coding_present']}",
        f"- all_expected_donor_specific_present: {summary['overall_flags']['all_expected_donor_specific_present']}",
        '',
        'Per-variant checks are in rna_variant_presence_checks.tsv.',
    ]
    pathlib.Path('rna_variant_presence_report.md').write_text('\n'.join(report_lines) + '\n')


if __name__ == '__main__':
    main()
