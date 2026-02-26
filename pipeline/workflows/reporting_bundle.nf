nextflow.enable.dsl = 2

include { WGS_ASSOCIATION } from './wgs_association'
include { RNA_INTEGRATION } from './rna_integration'
include { REPORT_BUNDLE } from '../modules/local/report_bundle/main'
include { REPORT_DASHBOARD } from '../modules/local/report_dashboard/main'

workflow REPORTING_BUNDLE {
  take:
    manifest_path

  main:
    assoc = WGS_ASSOCIATION(manifest_path)
    rna = RNA_INTEGRATION(manifest_path)
    report = REPORT_BUNDLE(
      assoc.assoc_diagnostics_summary,
      assoc.assoc_qc_summary,
      assoc.assoc_qc_top_hits,
      rna.rna_assoc_qc_summary,
      rna.rna_assoc_qc_top_hits,
      rna.rna_variant_presence_summary,
      rna.rna_ingest_report,
    )
    dashboard = REPORT_DASHBOARD(report.metrics_json, report.top_signals)

  emit:
    report_bundle_summary = report.summary_md
    report_bundle_metrics = report.metrics_json
    report_bundle_top_signals = report.top_signals
    report_dashboard_html = dashboard.html
    report_dashboard_data = dashboard.data
    assoc_qc_report = assoc.assoc_qc_report
    rna_assoc_qc_report = rna.rna_assoc_qc_report
}
