process EBV_VARIANT_MATRIX {
  label 'variant_compute'

  input:
    path sample_variant_tsvs

  output:
    path 'ebv_variant_catalog.tsv', emit: variant_catalog
    path 'ebv_variant_af_matrix.tsv', emit: af_matrix
    path 'ebv_variant_binary_matrix.tsv', emit: binary_matrix
    path 'ebv_variant_matrix_summary.json', emit: summary

  script:
  """
  python3 ${projectDir}/scripts/build_ebv_variant_matrices.py ${params.variant_binary_af_threshold}
  """
}
