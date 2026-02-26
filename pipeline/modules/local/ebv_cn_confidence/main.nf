process EBV_CN_CONFIDENCE {
  label 'cn_compute'

  input:
    path normalized_table
    path cohort_summary

  output:
    path 'ebv_cn_confidence.tsv', emit: confidence_table
    path 'ebv_cn_confidence_summary.json', emit: confidence_summary

  script:
  """
  python3 ${projectDir}/scripts/build_ebv_cn_confidence.py ${normalized_table} ${cohort_summary}
  """
}
