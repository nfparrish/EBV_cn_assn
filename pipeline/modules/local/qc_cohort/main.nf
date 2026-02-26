process QC_COHORT {
  label 'status_tracking'

  input:
    path qc_jsons

  output:
    path 'qc_cohort_metrics.tsv', emit: metrics_tsv
    path 'qc_cohort_summary.json', emit: summary_json
    path 'qc_cohort_report.md', emit: report_md

  script:
  """
  python3 ${projectDir}/scripts/build_qc_cohort.py
  """
}
