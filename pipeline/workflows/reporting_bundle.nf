nextflow.enable.dsl = 2

include { WGS_ASSOCIATION } from './wgs_association'
include { RNA_INTEGRATION } from './rna_integration'
include { REPORT_BUNDLE } from '../modules/local/report_bundle/main'
include { REPORT_DASHBOARD } from '../modules/local/report_dashboard/main'
include { REPORT_REPRO_CHECK } from '../modules/local/report_repro_check/main'

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
    repro = REPORT_REPRO_CHECK(
      report.summary_md,
      report.metrics_json,
      report.top_signals,
      dashboard.html,
      dashboard.data,
    )

  emit:
    report_bundle_summary = report.summary_md
    report_bundle_metrics = report.metrics_json
    report_bundle_top_signals = report.top_signals
    report_dashboard_html = dashboard.html
    report_dashboard_data = dashboard.data
    report_repro_manifest = repro.manifest
    report_repro_checksums = repro.checksums
    report_repro_report = repro.report
    assoc_qc_report = assoc.assoc_qc_report
    rna_assoc_qc_report = rna.rna_assoc_qc_report
}
