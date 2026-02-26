process REPORT_DASHBOARD {
  label 'assoc_compute'
  tag 'report_dashboard'

  input:
    path report_metrics
    path report_top_signals

  output:
    path 'report_dashboard.html', emit: html
    path 'report_dashboard_data.json', emit: data

  script:
  """
  python3 ${projectDir}/scripts/build_report_dashboard.py \
    ${report_metrics} \
    ${report_top_signals}
  """
}
