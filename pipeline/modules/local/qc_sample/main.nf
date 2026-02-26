process QC_SAMPLE {
  label 'staging'
  tag { sample_id }

  input:
    tuple val(sample_id), path(staged_json), path(status_json)

  output:
    tuple val(sample_id), path("${sample_id}.qc.json"), emit: qc

  script:
  """
  python3 - <<'PY'
  import json
  import pathlib
  import random

  sample_id = "${sample_id}"
  staged = json.loads(pathlib.Path("${staged_json}").read_text())
  status = json.loads(pathlib.Path("${status_json}").read_text())

  seed = sum(ord(c) for c in sample_id)
  rng = random.Random(seed)

  simulated_read_count = rng.randint(20_000_000, 60_000_000)
  simulated_q30_rate = round(rng.uniform(0.84, 0.97), 4)
  simulated_mapping_rate = round(rng.uniform(0.90, 0.995), 4)
  simulated_mean_coverage = round(rng.uniform(18.0, 42.0), 2)

  payload = {
      "sample_id": sample_id,
      "status": status.get("status"),
      "dataset": staged.get("dataset"),
      "population": staged.get("population"),
      "sequence_lab": staged.get("sequence_lab"),
      "metrics": {
          "read_count": simulated_read_count,
          "q30_rate": simulated_q30_rate,
          "mapping_rate": simulated_mapping_rate,
          "mean_coverage": simulated_mean_coverage,
      },
      "warnings": [
          "low_coverage" if simulated_mean_coverage < 20 else None,
          "low_q30" if simulated_q30_rate < 0.88 else None,
      ]
  }
  payload["warnings"] = [w for w in payload["warnings"] if w]

  pathlib.Path(f"{sample_id}.qc.json").write_text(json.dumps(payload, indent=2))
  PY
  """
}
