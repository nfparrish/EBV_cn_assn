process VALIDATE_MANIFEST {
  label 'manifest_validation'
  tag { manifest_file.baseName }

  input:
    path manifest_file

  output:
    path 'validated_manifest.csv', emit: validated_manifest
    path 'manifest_validation_report.json', emit: report

  script:
  """
  python3 - <<'PY'
  import csv
  import json
  import pathlib
  import sys

  manifest_path = pathlib.Path("${manifest_file}")
  required = ["sample_id", "dataset", "population", "sex", "sequence_lab"]
  allowed_dataset = {"1000G_PHASE3", "OTHER_PUBLIC"}
  allowed_sex = {"M", "F", "U"}

  with manifest_path.open("r", newline="") as handle:
      first_line = handle.readline()
      delimiter = '\t' if '\t' in first_line else ','

  with manifest_path.open("r", newline="") as handle:
      reader = csv.DictReader(handle, delimiter=delimiter)
      header = reader.fieldnames or []
      missing_columns = [field for field in required if field not in header]
      if missing_columns:
          raise SystemExit(f"Missing required manifest columns: {', '.join(missing_columns)}")

      rows = []
      errors = []
      seen = set()

      for row_idx, row in enumerate(reader, start=2):
          normalized = {k: (v or '').strip() for k, v in row.items()}
          sample_id = normalized.get("sample_id", "")

          if not sample_id:
              errors.append({"row": row_idx, "error": "sample_id is required"})
              continue

          if sample_id in seen:
              errors.append({"row": row_idx, "sample_id": sample_id, "error": "duplicate sample_id"})
              continue
          seen.add(sample_id)

          if normalized["dataset"] not in allowed_dataset:
              errors.append({"row": row_idx, "sample_id": sample_id, "error": "dataset must be 1000G_PHASE3 or OTHER_PUBLIC"})
              continue

          if normalized["sex"] not in allowed_sex:
              errors.append({"row": row_idx, "sample_id": sample_id, "error": "sex must be one of M/F/U"})
              continue

          has_pair = bool(normalized.get("read1_path")) and bool(normalized.get("read2_path"))
          has_bam = bool(normalized.get("bam_path"))
          has_cram = bool(normalized.get("cram_path"))
          if not (has_pair or has_bam or has_cram):
              errors.append({"row": row_idx, "sample_id": sample_id, "error": "must provide read1+read2 or bam_path or cram_path"})
              continue

          rows.append(normalized)

  if errors:
      report = {
          "manifest": str(manifest_path),
          "status": "failed",
          "valid_rows": len(rows),
          "invalid_rows": len(errors),
          "errors": errors,
      }
      pathlib.Path("manifest_validation_report.json").write_text(json.dumps(report, indent=2))
      raise SystemExit("Manifest validation failed; see manifest_validation_report.json")

  output_cols = list(header)
  for row in rows:
      for key in row.keys():
          if key not in output_cols:
              output_cols.append(key)
  with pathlib.Path("validated_manifest.csv").open("w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=output_cols)
      writer.writeheader()
      writer.writerows(rows)

  report = {
      "manifest": str(manifest_path),
      "status": "passed",
      "rows": len(rows),
      "columns": output_cols,
  }
  pathlib.Path("manifest_validation_report.json").write_text(json.dumps(report, indent=2))
  PY
  """
}
