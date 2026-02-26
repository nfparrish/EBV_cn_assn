process STAGE_INPUTS {
  label 'staging'
  tag { row.sample_id }

  input:
    val row

  output:
    tuple val(row.sample_id), path("${row.sample_id}.inputs.json"), emit: staged

  script:
  """
  python3 - <<'PY'
  import json
  import pathlib

  row = json.loads('''${groovy.json.JsonOutput.toJson(row)}''')
  sample_id = row["sample_id"]

  payload = {
      "sample_id": sample_id,
      "dataset": row.get("dataset"),
      "population": row.get("population"),
      "sex": row.get("sex"),
      "sequence_lab": row.get("sequence_lab"),
      "read1_path": row.get("read1_path"),
      "read2_path": row.get("read2_path"),
      "bam_path": row.get("bam_path"),
      "cram_path": row.get("cram_path"),
      "md5": row.get("md5"),
      "source_exists": {
          "read1_path": bool(row.get("read1_path")) and pathlib.Path(row.get("read1_path")).exists(),
          "read2_path": bool(row.get("read2_path")) and pathlib.Path(row.get("read2_path")).exists(),
          "bam_path": bool(row.get("bam_path")) and pathlib.Path(row.get("bam_path")).exists(),
          "cram_path": bool(row.get("cram_path")) and pathlib.Path(row.get("cram_path")).exists(),
      }
  }

  pathlib.Path(f"{sample_id}.inputs.json").write_text(json.dumps(payload, indent=2))
  PY
  """
}
