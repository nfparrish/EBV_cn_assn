process VARIANT_PROVENANCE {
  label 'variant_compute'
  tag 'variant_provenance'

  input:
    path filter_counts
    path filter_rules
    path ebv_matrix_summary
    path host_ingest_report
    path host_harmonize_report

  output:
    path 'variant_provenance_summary.json', emit: summary
    path 'variant_provenance_counts.tsv', emit: counts
    path 'variant_provenance_report.md', emit: report

  script:
  """
  python3 ${projectDir}/scripts/build_variant_provenance.py \
    ${filter_counts} \
    ${filter_rules} \
    ${ebv_matrix_summary} \
    ${host_ingest_report} \
    ${host_harmonize_report}
  """
}
