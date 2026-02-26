process REPORT_REPRO_CHECK {
  label 'assoc_compute'
  tag 'report_repro_check'

  input:
    path report_bundle_summary
    path report_bundle_metrics
    path report_bundle_top_signals
    path report_dashboard_html
    path report_dashboard_data

  output:
    path 'report_reproducibility_manifest.json', emit: manifest
    path 'report_reproducibility_checksums.tsv', emit: checksums
    path 'report_reproducibility_report.md', emit: report

  script:
  """
  python3 ${projectDir}/scripts/build_report_reproducibility.py \
    ${report_bundle_summary} \
    ${report_bundle_metrics} \
    ${report_bundle_top_signals} \
    ${report_dashboard_html} \
    ${report_dashboard_data}
  """
}
