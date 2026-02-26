process ASSOC_QC_SUMMARY {
  label 'assoc_compute'
  tag 'assoc_qc_summary'

  input:
    path diagnostics_summary
    path ebv_results
    path host_results

  output:
    path 'association_qc_summary.json', emit: qc_summary
    path 'association_qc_report.md', emit: qc_report
    path 'association_qc_top_hits.tsv', emit: top_hits

  script:
  """
  python3 ${projectDir}/scripts/build_assoc_qc_summary.py ${diagnostics_summary} ${ebv_results} ${host_results}
  """
}
