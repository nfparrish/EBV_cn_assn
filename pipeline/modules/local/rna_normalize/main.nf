process RNA_NORMALIZE {
  label 'rna_compute'
  tag 'rna_normalize'

  input:
    path expression_canonical

  output:
    path 'rna_expression_log2tpm.tsv', emit: expression_log2
    path 'rna_expression_quantile.tsv', emit: expression_quantile
    path 'rna_normalization_report.json', emit: normalization_report

  script:
  """
  python3 ${projectDir}/scripts/normalize_rna_expression.py ${expression_canonical}
  """
}
