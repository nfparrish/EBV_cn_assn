import csv
import json
import math
import pathlib
import sys

import numpy as np


def read_tsv(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader)


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
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


def build_covariates(rows):
    sample_ids = [row['sample_id'] for row in rows]

    pcs = []
    for row in rows:
        pcs.append([
            to_float(row.get('pc1')),
            to_float(row.get('pc2')),
            to_float(row.get('pc3')),
            to_float(row.get('pc4')),
            to_float(row.get('pc5')),
        ])

    sex_vec = [1.0 if (row.get('sex', '').upper() == 'F') else 0.0 for row in rows]

    sequence_labs = [row.get('sequence_lab', '') for row in rows]
    batches = [row.get('batch', '') for row in rows]
    lab_batches = [row.get('lab_batch', '') for row in rows]

    has_separate = any(sequence_labs) and any(batches)

    cov_matrix = []
    headers = ['intercept', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'sex_F']

    unique_labs = sorted({x for x in sequence_labs if x})
    unique_batches = sorted({x for x in batches if x})
    unique_lab_batches = sorted({x for x in lab_batches if x})

    lab_cols = unique_labs[1:] if len(unique_labs) > 1 else []
    batch_cols = unique_batches[1:] if len(unique_batches) > 1 else []
    lab_batch_cols = unique_lab_batches[1:] if len(unique_lab_batches) > 1 else []

    if has_separate:
        headers += [f'lab::{x}' for x in lab_cols]
        headers += [f'batch::{x}' for x in batch_cols]
    else:
        headers += [f'lab_batch::{x}' for x in lab_batch_cols]

    for idx, sid in enumerate(sample_ids):
        row = [1.0] + pcs[idx] + [sex_vec[idx]]

        if has_separate:
            row += [1.0 if sequence_labs[idx] == lab else 0.0 for lab in lab_cols]
            row += [1.0 if batches[idx] == batch else 0.0 for batch in batch_cols]
        else:
            row += [1.0 if lab_batches[idx] == lb else 0.0 for lb in lab_batch_cols]

        cov_matrix.append(row)

    return np.asarray(cov_matrix, dtype=float), headers, ('separate' if has_separate else 'combined')


def residualize(vector, covariates):
    pinv = np.linalg.pinv(covariates)
    fitted = covariates @ (pinv @ vector)
    return vector - fitted


def run_variant_assoc(pathway, assoc_rows, variant_rows, min_samples):
    y = np.asarray([to_float(row.get('cn_zscore')) for row in assoc_rows], dtype=float)
    cov_matrix, cov_headers, cov_mode = build_covariates(assoc_rows)

    sample_ids = [row['sample_id'] for row in assoc_rows]

    by_variant = {}
    for row in variant_rows:
        sid = (row.get('sample_id') or '').strip()
        if not sid:
            continue

        if pathway == 'ebv':
            variant_id = f"{row.get('chrom')}:{row.get('pos')}:{row.get('ref')}>{row.get('alt')}"
            value = to_float(row.get('af'))
        else:
            variant_id = row.get('variant_id', '').strip()
            value = to_float(row.get('dosage'))

        if not variant_id:
            continue

        by_variant.setdefault(variant_id, {})[sid] = value

    results = []

    for variant_id in sorted(by_variant.keys()):
        sample_to_val = by_variant[variant_id]
        x_raw = np.asarray([sample_to_val.get(sid, 0.0) for sid in sample_ids], dtype=float)

        non_zero = int(np.count_nonzero(x_raw))
        if len(sample_ids) < min_samples or np.var(x_raw) <= 1e-12:
            results.append({
                'variant_id': variant_id,
                'pathway': pathway,
                'n_samples': len(sample_ids),
                'n_nonzero_predictor': non_zero,
                'effect_size': 'nan',
                'std_error': 'nan',
                'z_score': 'nan',
                'p_value': 'nan',
                'p_value_fdr': 'nan',
                'model_family': 'gaussian',
                'model_link': 'identity',
                'covariate_mode': cov_mode,
            })
            continue

        y_res = residualize(y, cov_matrix)
        x_res = residualize(x_raw, cov_matrix)

        denom = float(np.dot(x_res, x_res))
        if denom <= 1e-12:
            results.append({
                'variant_id': variant_id,
                'pathway': pathway,
                'n_samples': len(sample_ids),
                'n_nonzero_predictor': non_zero,
                'effect_size': 'nan',
                'std_error': 'nan',
                'z_score': 'nan',
                'p_value': 'nan',
                'p_value_fdr': 'nan',
                'model_family': 'gaussian',
                'model_link': 'identity',
                'covariate_mode': cov_mode,
            })
            continue

        beta = float(np.dot(x_res, y_res) / denom)
        residual = y_res - beta * x_res

        dof = max(len(sample_ids) - cov_matrix.shape[1] - 1, 1)
        sigma2 = float(np.dot(residual, residual) / dof)
        std_error = math.sqrt(max(sigma2 / denom, 0.0))

        if std_error > 0:
            z_score = beta / std_error
            p_value = math.erfc(abs(z_score) / math.sqrt(2.0))
        else:
            z_score = 0.0
            p_value = 1.0

        results.append({
            'variant_id': variant_id,
            'pathway': pathway,
            'n_samples': len(sample_ids),
            'n_nonzero_predictor': non_zero,
            'effect_size': f"{beta:.8g}",
            'std_error': f"{std_error:.8g}",
            'z_score': f"{z_score:.8g}",
            'p_value': f"{p_value:.8g}",
            'p_value_fdr': None,
            'model_family': 'gaussian',
            'model_link': 'identity',
            'covariate_mode': cov_mode,
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

    return results, cov_headers, cov_mode


def main():
    assoc_path = pathlib.Path(sys.argv[1])
    variant_path = pathlib.Path(sys.argv[2])
    pathway = sys.argv[3]
    out_tsv = pathlib.Path(sys.argv[4])
    out_json = pathlib.Path(sys.argv[5])
    min_samples = int(sys.argv[6]) if len(sys.argv) > 6 else 3

    assoc_rows = read_tsv(assoc_path)
    variant_rows = read_tsv(variant_path)

    results, cov_headers, cov_mode = run_variant_assoc(pathway, assoc_rows, variant_rows, min_samples)

    header = [
        'variant_id', 'pathway', 'n_samples', 'n_nonzero_predictor',
        'effect_size', 'std_error', 'z_score', 'p_value', 'p_value_fdr',
        'model_family', 'model_link', 'covariate_mode',
    ]

    with out_tsv.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    finite_p = []
    for row in results:
        try:
            value = float(row['p_value'])
            if not math.isnan(value) and not math.isinf(value):
                finite_p.append(value)
        except Exception:
            pass

    report = {
        'pathway': pathway,
        'result_count': len(results),
        'min_samples_required': min_samples,
        'covariate_mode': cov_mode,
        'covariate_terms': cov_headers,
        'min_p_value': min(finite_p) if finite_p else None,
    }
    out_json.write_text(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
