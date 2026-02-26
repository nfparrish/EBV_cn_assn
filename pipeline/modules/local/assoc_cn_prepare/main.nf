process ASSOC_CN_PREPARE {
  label 'assoc_compute'
  tag 'assoc_cn_prepare'

  input:
    path cn_table
    path covariate_table

  output:
    path 'association_input.tsv', emit: model_input
    path 'association_cn_prepare_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/prepare_assoc_cn.py ${cn_table} ${covariate_table} ${params.association_cn_column}
  """
}
