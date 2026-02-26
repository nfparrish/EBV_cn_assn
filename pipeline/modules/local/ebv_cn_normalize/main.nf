process EBV_CN_NORMALIZE {
  label 'cn_compute'

  input:
    path cn_raw_jsons

  output:
    path 'ebv_cn_raw_cohort.tsv', emit: raw_table
    path 'ebv_cn_normalized.tsv', emit: normalized_table
    path 'ebv_cn_cohort_summary.json', emit: summary

  script:
  """
  python3 ${projectDir}/scripts/normalize_ebv_cn.py
  """
}
