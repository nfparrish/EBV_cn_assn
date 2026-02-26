process REPORT_BUNDLE {
  label 'assoc_compute'
  tag 'report_bundle'

  input:
    path assoc_diagnostics_summary
    path assoc_qc_summary
    path assoc_qc_top_hits
    path rna_assoc_qc_summary
    path rna_assoc_qc_top_hits
    path rna_variant_presence_summary
    path rna_ingest_report

  output:
    path 'report_bundle_summary.md', emit: summary_md
    path 'report_bundle_metrics.json', emit: metrics_json
    path 'report_bundle_top_signals.tsv', emit: top_signals

  script:
  """
  python3 ${projectDir}/scripts/build_report_bundle.py \
    ${assoc_diagnostics_summary} \
    ${assoc_qc_summary} \
    ${assoc_qc_top_hits} \
    ${rna_assoc_qc_summary} \
    ${rna_assoc_qc_top_hits} \
    ${rna_variant_presence_summary} \
    ${rna_ingest_report}
  """
}
