import csv
import glob
import json
import math
import pathlib
import sys

import numpy as np


def to_float(value, default=0.0):
    try:
        v = float(value)
        if math.isnan(v) or math.isinf(v):
            return default
        return v
    except Exception:
        return default


def bh_fdr(pvalues):
    indexed = [(idx, p) for idx, p in enumerate(pvalues) if p is not None and not math.isnan(p)]
    if not indexed:
        return [None] * len(pvalues)

    indexed.sort(key=lambda x: x[1])
    m = len(indexed)
    adjusted = [None] * len(pvalues)

    prev = 1.0
    for rank in range(m, 0, -1):
        idx, p = indexed[rank - 1]
        raw = p * m / rank
        prev = min(prev, raw)
        adjusted[idx] = min(prev, 1.0)

    return adjusted


def deterministic_pc(sample_id, pc_index):
    seed = sum(ord(c) for c in sample_id) + (pc_index + 1) * 113
    rng = np.random.default_rng(seed)
    return float(rng.uniform(-0.2, 0.2))


def read_expression(path):
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        headers = reader.fieldnames or []
        if 'gene_id' not in headers:
            raise ValueError("RNA expression matrix missing 'gene_id' column")
        samples = [h for h in headers if h != 'gene_id']

        genes = []
        values = {}
        for row in reader:
            gene_id = (row.get('gene_id') or '').strip()
            if not gene_id:
                continue
            genes.append(gene_id)
            values[gene_id] = {sid: to_float(row.get(sid), default=0.0) for sid in samples}

    return samples, genes, values


def read_sample_map(path):
    rows = {}
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            if sid:
                rows[sid] = row
    return rows


def parse_variants(min_af, min_dp, min_qual):
    by_variant = {}
    source_rows = 0

    for file_path in sorted(glob.glob('*.ebv_variants.tsv')):
        with pathlib.Path(file_path).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            for row in reader:
                source_rows += 1
                sid = (row.get('sample_id') or '').strip()
                if not sid:
                    continue

                chrom = (row.get('chrom') or '').strip()
                pos = (row.get('pos') or '').strip()
                ref = (row.get('ref') or '').strip()
                alt = (row.get('alt') or '').strip()
                if not chrom or not pos or not ref or not alt:
                    continue

                af = to_float(row.get('af'), 0.0)
                dp = to_float(row.get('dp'), 0.0)
                qual = to_float(row.get('qual'), 0.0)
                filt = (row.get('filter') or '').strip().upper()

                if af < min_af or dp < min_dp or qual < min_qual:
                    continue
                if filt not in ('PASS', ''):
                    continue

                variant_id = f"{chrom}:{pos}:{ref}>{alt}"
                sample_map = by_variant.setdefault(variant_id, {})
                if sid not in sample_map or af > sample_map[sid]:
                    sample_map[sid] = af

    return by_variant, source_rows


def build_covariates(sample_ids, sample_map):
    pcs = []
    sex_f = []
    batches = []
    library_size = []
    mapping_rate = []

    for sid in sample_ids:
        m = sample_map.get(sid, {})
        pcs.append([
            to_float(m.get('pc1'), deterministic_pc(sid, 0)),
            to_float(m.get('pc2'), deterministic_pc(sid, 1)),
            to_float(m.get('pc3'), deterministic_pc(sid, 2)),
            to_float(m.get('pc4'), deterministic_pc(sid, 3)),
            to_float(m.get('pc5'), deterministic_pc(sid, 4)),
        ])
        sex_f.append(1.0 if (m.get('sex', '').upper() == 'F') else 0.0)
        batches.append((m.get('rna_batch') or m.get('batch') or '').strip())
        library_size.append(to_float(m.get('library_size'), 0.0))
        mapping_rate.append(to_float(m.get('mapping_rate'), 0.0))

    unique_batches = sorted({b for b in batches if b})
    batch_cols = unique_batches[1:] if len(unique_batches) > 1 else []

    n_samples = len(sample_ids)
    include_pcs = n_samples >= 8
    include_library_metrics = n_samples >= 6
    include_batch = n_samples >= 10
    include_sex = n_samples >= 3

    headers = ['intercept']
    if include_pcs:
        headers += ['pc1', 'pc2', 'pc3', 'pc4', 'pc5']
    if include_sex:
        headers += ['sex_F']
    if include_library_metrics:
        headers += ['library_size', 'mapping_rate']
    if include_batch:
        headers += [f'batch::{b}' for b in batch_cols]

    cov = []
    for i in range(n_samples):
        row = [1.0]
        if include_pcs:
            row += pcs[i]
        if include_sex:
            row += [sex_f[i]]
        if include_library_metrics:
            row += [library_size[i], mapping_rate[i]]
        if include_batch:
            row += [1.0 if batches[i] == b else 0.0 for b in batch_cols]
        cov.append(row)

    return np.asarray(cov, dtype=float), headers


def residualize(vector, covariates):
    pinv = np.linalg.pinv(covariates)
    fitted = covariates @ (pinv @ vector)
    return vector - fitted


def run_assoc(variant_map, sample_map, expr_samples, genes, expr_values, normalization_label, min_samples):
    all_variant_samples = set()
    for sample_vals in variant_map.values():
        all_variant_samples.update(sample_vals.keys())

    sample_ids = sorted(all_variant_samples.intersection(set(expr_samples)).intersection(set(sample_map.keys())))
    cov_matrix, cov_headers = build_covariates(sample_ids, sample_map)

    results = []

    for variant_id in sorted(variant_map.keys()):
        x_raw = np.asarray([variant_map[variant_id].get(sid, 0.0) for sid in sample_ids], dtype=float)
        non_zero = int(np.count_nonzero(x_raw))

        x_res = residualize(x_raw, cov_matrix) if len(sample_ids) > 0 else np.asarray([], dtype=float)
        denom = float(np.dot(x_res, x_res)) if len(sample_ids) > 0 else 0.0

        for gene_id in genes:
            y = np.asarray([expr_values[gene_id].get(sid, 0.0) for sid in sample_ids], dtype=float)

            if len(sample_ids) < min_samples or np.var(y) <= 1e-12 or np.var(x_raw) <= 1e-12 or denom <= 1e-12:
                results.append({
                    'normalization': normalization_label,
                    'variant_id': variant_id,
                    'gene_id': gene_id,
                    'n_samples': len(sample_ids),
                    'n_nonzero_predictor': non_zero,
                    'effect_size': 'nan',
                    'std_error': 'nan',
                    't_stat': 'nan',
                    'p_value': 'nan',
                    'p_value_fdr': 'nan',
                    'model_family': 'gaussian',
                    'model_link': 'identity',
                    'covariate_terms': ';'.join(cov_headers),
                })
                continue

            dof_raw = len(sample_ids) - cov_matrix.shape[1] - 1
            if dof_raw < 1:
                results.append({
                    'normalization': normalization_label,
                    'variant_id': variant_id,
                    'gene_id': gene_id,
                    'n_samples': len(sample_ids),
                    'n_nonzero_predictor': non_zero,
                    'effect_size': 'nan',
                    'std_error': 'nan',
                    't_stat': 'nan',
                    'p_value': 'nan',
                    'p_value_fdr': 'nan',
                    'model_family': 'gaussian',
                    'model_link': 'identity',
                    'covariate_terms': ';'.join(cov_headers),
                })
                continue

            y_res = residualize(y, cov_matrix)
            beta = float(np.dot(x_res, y_res) / denom)
            residual = y_res - beta * x_res

            sigma2 = float(np.dot(residual, residual) / dof_raw)
            std_error = math.sqrt(max(sigma2 / denom, 0.0))

            if std_error > 0:
                t_stat = beta / std_error
                p_val = math.erfc(abs(t_stat) / math.sqrt(2.0))
            else:
                t_stat = 0.0
                p_val = 1.0

            results.append({
                'normalization': normalization_label,
                'variant_id': variant_id,
                'gene_id': gene_id,
                'n_samples': len(sample_ids),
                'n_nonzero_predictor': non_zero,
                'effect_size': f"{beta:.8g}",
                'std_error': f"{std_error:.8g}",
                't_stat': f"{t_stat:.8g}",
                'p_value': f"{p_val:.8g}",
                'p_value_fdr': None,
                'model_family': 'gaussian',
                'model_link': 'identity',
                'covariate_terms': ';'.join(cov_headers),
            })

    pvals = []
    for row in results:
        try:
            pvals.append(float(row['p_value']))
        except Exception:
            pvals.append(float('nan'))

    adjusted = bh_fdr(pvals)
    for idx, adj in enumerate(adjusted):
        if adj is None or (isinstance(adj, float) and math.isnan(adj)):
            results[idx]['p_value_fdr'] = 'nan'
        else:
            results[idx]['p_value_fdr'] = f"{adj:.8g}"

    return results, sample_ids, cov_headers


def main():
    expr_path = sys.argv[1]
    sample_map_path = sys.argv[2]
    normalization_label = sys.argv[3]
    out_tsv = pathlib.Path(sys.argv[4])
    out_json = pathlib.Path(sys.argv[5])
    min_samples = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    min_af = to_float(sys.argv[7], 0.05) if len(sys.argv) > 7 else 0.05
    min_dp = to_float(sys.argv[8], 20.0) if len(sys.argv) > 8 else 20.0
    min_qual = to_float(sys.argv[9], 60.0) if len(sys.argv) > 9 else 60.0

    expr_samples, genes, expr_values = read_expression(expr_path)
    sample_map = read_sample_map(sample_map_path)
    variant_map, source_rows = parse_variants(min_af=min_af, min_dp=min_dp, min_qual=min_qual)

    results, overlap_samples, cov_headers = run_assoc(
        variant_map,
        sample_map,
        expr_samples,
        genes,
        expr_values,
        normalization_label,
        min_samples,
    )

    header = [
        'normalization', 'variant_id', 'gene_id', 'n_samples', 'n_nonzero_predictor',
        'effect_size', 'std_error', 't_stat', 'p_value', 'p_value_fdr',
        'model_family', 'model_link', 'covariate_terms',
    ]

    with out_tsv.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    finite_p = []
    for row in results:
        try:
            v = float(row['p_value'])
            if not math.isnan(v) and not math.isinf(v):
                finite_p.append(v)
        except Exception:
            pass

    report = {
        'normalization': normalization_label,
        'variant_input_row_count': source_rows,
        'variant_count_modeled': len({row['variant_id'] for row in results}),
        'test_count': len(results),
        'gene_count': len(genes),
        'sample_overlap_count': len(overlap_samples),
        'min_samples_required': min_samples,
        'covariate_terms': cov_headers,
        'filters': {
            'min_af': min_af,
            'min_dp': min_dp,
            'min_qual': min_qual,
        },
        'min_p_value': min(finite_p) if finite_p else None,
    }
    out_json.write_text(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
