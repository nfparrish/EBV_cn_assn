import csv
import json
import math
import pathlib
import sys


def read_matrix(path):
    with pathlib.Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        headers = reader.fieldnames or []
        if 'gene_id' not in headers:
            raise ValueError("Canonical matrix missing 'gene_id'")
        samples = [h for h in headers if h != 'gene_id']

        genes = []
        matrix = []
        for row in reader:
            genes.append(row['gene_id'])
            matrix.append([float(row.get(s, 0.0) or 0.0) for s in samples])

    return genes, samples, matrix


def write_matrix(path, genes, samples, matrix):
    with pathlib.Path(path).open('w') as handle:
        handle.write('\t'.join(['gene_id'] + samples) + '\n')
        for gene, values in zip(genes, matrix):
            handle.write('\t'.join([gene] + [f"{v:.8f}" for v in values]) + '\n')


def log2_tpm_matrix(matrix):
    out = []
    for row in matrix:
        out.append([math.log2(v + 1.0) for v in row])
    return out


def quantile_normalize(matrix):
    n_genes = len(matrix)
    if n_genes == 0:
        return []
    n_samples = len(matrix[0])

    columns = [[] for _ in range(n_samples)]
    for row in matrix:
        for j, v in enumerate(row):
            columns[j].append(v)

    sorted_cols = [sorted(col) for col in columns]
    rank_means = [sum(sorted_cols[j][i] for j in range(n_samples)) / n_samples for i in range(n_genes)]

    normalized_cols = []
    for col in columns:
        order = sorted(range(n_genes), key=lambda i: col[i])
        out_col = [0.0] * n_genes
        for rank, gene_idx in enumerate(order):
            out_col[gene_idx] = rank_means[rank]
        normalized_cols.append(out_col)

    normalized = [[0.0] * n_samples for _ in range(n_genes)]
    for j in range(n_samples):
        for i in range(n_genes):
            normalized[i][j] = normalized_cols[j][i]

    return normalized


canonical_path = sys.argv[1]

genes, samples, matrix = read_matrix(canonical_path)

log2_matrix = log2_tpm_matrix(matrix)
quant_matrix = quantile_normalize(matrix)

write_matrix('rna_expression_log2tpm.tsv', genes, samples, log2_matrix)
write_matrix('rna_expression_quantile.tsv', genes, samples, quant_matrix)

report = {
    'sample_count': len(samples),
    'gene_count': len(genes),
    'methods': [
        {'name': 'log2_tpm_plus_1', 'description': 'log2(TPM + 1) transformation'},
        {'name': 'quantile_normalized', 'description': 'sample-wise quantile normalization'},
    ],
    'matrix_dimensions': {'rows_genes': len(genes), 'cols_samples': len(samples)},
}
pathlib.Path('rna_normalization_report.json').write_text(json.dumps(report, indent=2))
