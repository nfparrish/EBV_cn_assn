process SAMPLE_STATUS {
  label 'status_tracking'
  tag { sample_id }

  input:
    tuple val(sample_id), path(staged_json)

  output:
    tuple val(sample_id), path("${sample_id}.status.json"), emit: status

  script:
  """
  python3 - <<'PY'
  import csv
  import json
  import os
  import pathlib

  sample_id = "${sample_id}"
  staged_path = pathlib.Path("${staged_json}")
  staged_payload = json.loads(staged_path.read_text())

  registry = "${params.status_registry ?: ''}".strip()
  prior_status = None

  if registry and pathlib.Path(registry).exists():
      with pathlib.Path(registry).open("r", newline="") as handle:
          reader = csv.DictReader(handle)
          for row in reader:
              if row.get("sample_id") == sample_id:
                  prior_status = row.get("status")
                  break

  if prior_status == "processed":
      status = "skipped"
      reason = "already_processed"
  else:
      status = "processed"
      reason = "new_or_retry"

  payload = {
      "sample_id": sample_id,
      "status": status,
      "reason": reason,
      "dataset": staged_payload.get("dataset"),
      "population": staged_payload.get("population"),
      "sequence_lab": staged_payload.get("sequence_lab"),
  }

  pathlib.Path(f"{sample_id}.status.json").write_text(json.dumps(payload, indent=2))
  PY
  """

}
