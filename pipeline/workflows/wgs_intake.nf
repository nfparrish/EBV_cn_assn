nextflow.enable.dsl = 2

include { VALIDATE_MANIFEST } from '../modules/local/validate_manifest/main'
include { STAGE_INPUTS } from '../modules/local/stage_inputs/main'
include { SAMPLE_STATUS } from '../modules/local/sample_status/main'
include { QC_SAMPLE } from '../modules/local/qc_sample/main'
include { QC_COHORT } from '../modules/local/qc_cohort/main'

process BUILD_STATUS_SUMMARY {
  label 'status_tracking'

  input:
    path status_jsons

  output:
    path 'sample_status_registry.tsv', emit: registry
    path 'sample_status_summary.json', emit: summary

  script:
  """
  python3 ${projectDir}/scripts/build_status_summary.py
  """
}

workflow WGS_INTAKE {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)

    validation = VALIDATE_MANIFEST(manifest_file)

    manifest_rows = validation.validated_manifest.splitCsv(header: true)
    staged = STAGE_INPUTS(manifest_rows)
    statuses = SAMPLE_STATUS(staged.staged)
    qc_sample = QC_SAMPLE(staged.staged.join(statuses.status, by: 0))
    qc_cohort = QC_COHORT(qc_sample.qc.map { it[1] }.collect())

    BUILD_STATUS_SUMMARY(statuses.status.map { it[1] }.collect())

  emit:
    validated_manifest = validation.validated_manifest
    validation_report = validation.report
    status_registry = BUILD_STATUS_SUMMARY.out.registry
    status_summary = BUILD_STATUS_SUMMARY.out.summary
    qc_metrics = qc_cohort.metrics_tsv
    qc_summary = qc_cohort.summary_json
    qc_report = qc_cohort.report_md
}
