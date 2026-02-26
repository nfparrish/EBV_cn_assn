process EBV_CN_QC {
  label 'cn_compute'

  input:
    path confidence_table

  output:
    path 'ebv_cn_qc_summary.json', emit: qc_summary
    path 'ebv_cn_qc_report.md', emit: qc_report
    path 'ebv_cn_low_confidence.tsv', emit: low_confidence

  script:
  """
  python3 ${projectDir}/scripts/build_ebv_cn_qc.py ${confidence_table}
  """
}
