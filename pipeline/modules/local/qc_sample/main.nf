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
    import shutil
    import subprocess

  sample_id = "${sample_id}"
  staged = json.loads(pathlib.Path("${staged_json}").read_text())
  status = json.loads(pathlib.Path("${status_json}").read_text())

    def run_command(cmd):
      completed = subprocess.run(cmd, capture_output=True, text=True)
      if completed.returncode != 0:
        raise RuntimeError(
          f"Command failed: {' '.join(cmd)}\n"
          f"STDOUT: {completed.stdout.strip()}\n"
          f"STDERR: {completed.stderr.strip()}"
        )
      return completed.stdout

    read_count = None
    q30_rate = None
    mapping_rate = None
    mean_coverage = None

    alignment_path = (staged.get('cram_path') or '').strip() or (staged.get('bam_path') or '').strip()
    if alignment_path and pathlib.Path(alignment_path).exists() and shutil.which('samtools'):
      try:
        total = int(run_command(['samtools', 'view', '-c', alignment_path]).strip() or '0')
        idxstats = run_command(['samtools', 'idxstats', alignment_path]).splitlines()
        mapped = 0
        for line in idxstats:
          parts = line.split('\t')
          if len(parts) == 4 and parts[0] != '*':
            try:
              mapped += int(parts[2])
            except ValueError:
              pass
        read_count = total
        mapping_rate = (float(mapped) / float(total)) if total > 0 else None

        total_bases = 0
        weighted_depth = 0.0
        coverage = run_command(['samtools', 'coverage', alignment_path]).splitlines()
        for line in coverage:
          if not line or line.startswith('#'):
            continue
          parts = line.split('\t')
          if len(parts) < 7:
            continue
          try:
            start = int(parts[1])
            end = int(parts[2])
            depth = float(parts[6])
          except ValueError:
            continue
          length = max(end - start + 1, 0)
          total_bases += length
          weighted_depth += depth * length
        mean_coverage = (weighted_depth / total_bases) if total_bases > 0 else None
      except Exception:
        pass

  payload = {
      "sample_id": sample_id,
      "status": status.get("status"),
      "dataset": staged.get("dataset"),
      "population": staged.get("population"),
      "sequence_lab": staged.get("sequence_lab"),
      "metrics": {
          "read_count": read_count,
          "q30_rate": q30_rate,
          "mapping_rate": round(mapping_rate, 6) if mapping_rate is not None else None,
          "mean_coverage": round(mean_coverage, 6) if mean_coverage is not None else None,
      },
      "warnings": [
          "missing_alignment_metrics" if read_count is None else None,
          "low_coverage" if (mean_coverage is not None and mean_coverage < 20) else None,
          "low_mapping_rate" if (mapping_rate is not None and mapping_rate < 0.9) else None,
      ]
  }
  payload["warnings"] = [w for w in payload["warnings"] if w]

  pathlib.Path(f"{sample_id}.qc.json").write_text(json.dumps(payload, indent=2))
  PY
  """
}
