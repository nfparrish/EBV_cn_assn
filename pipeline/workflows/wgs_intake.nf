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
  python3 - <<'PY'
  import glob
  import json
  import pathlib

  rows = []
  for file_path in glob.glob('*.status.json'):
      rows.append(json.loads(pathlib.Path(file_path).read_text()))

  rows = sorted(rows, key=lambda item: item["sample_id"])

  header = ["sample_id", "status", "reason", "dataset", "population", "sequence_lab"]
  with pathlib.Path("sample_status_registry.tsv").open("w") as handle:
      handle.write("\t".join(header) + "\n")
      for row in rows:
          handle.write("\t".join(str(row.get(col, "")) for col in header) + "\n")

  summary = {
      "total": len(rows),
      "processed": sum(1 for row in rows if row.get("status") == "processed"),
      "skipped": sum(1 for row in rows if row.get("status") == "skipped"),
  }
  pathlib.Path("sample_status_summary.json").write_text(json.dumps(summary, indent=2))
  PY
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
    qc_cohort = QC_COHORT(qc_sample.qc.collect())

    BUILD_STATUS_SUMMARY(statuses.status.collect())

  emit:
    validated_manifest = validation.validated_manifest
    validation_report = validation.report
    status_registry = BUILD_STATUS_SUMMARY.out.registry
    status_summary = BUILD_STATUS_SUMMARY.out.summary
    qc_metrics = qc_cohort.metrics_tsv
    qc_summary = qc_cohort.summary_json
    qc_report = qc_cohort.report_md
}
