import json
import pathlib
import shutil
import subprocess
import sys


def run_command(cmd):
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDOUT: {completed.stdout.strip()}\n"
            f"STDERR: {completed.stderr.strip()}"
        )
    return completed.stdout


sample_id = sys.argv[1]
staged_path = pathlib.Path(sys.argv[2])
status_path = pathlib.Path(sys.argv[3])
cram_reference = sys.argv[4].strip() if len(sys.argv) > 4 else ''

staged = json.loads(staged_path.read_text())
status = json.loads(status_path.read_text())

read_count = None
q30_rate = None
mapping_rate = None
mean_coverage = None

alignment_path = (staged.get('cram_path') or '').strip() or (staged.get('bam_path') or '').strip()
if alignment_path and pathlib.Path(alignment_path).exists() and shutil.which('samtools'):
    try:
        ref_env = {'CRAM_REFERENCE': cram_reference} if cram_reference else {}
        import os
        env = os.environ.copy()
        env.update(ref_env)

        # Total reads via flagstat (inherits CRAM_REFERENCE from env)
        flagstat = subprocess.run(
            ['samtools', 'flagstat', '-@4', alignment_path],
            capture_output=True, text=True, env=env
        )
        total = 0
        mapped = 0
        if flagstat.returncode == 0:
            for line in flagstat.stdout.splitlines():
                if 'in total' in line:
                    total = int(line.split()[0])
                if 'mapped (' in line and 'primary' not in line:
                    mapped = int(line.split()[0])
        read_count = total if total > 0 else None
        mapping_rate = (float(mapped) / float(total)) if total > 0 else None

        # Mean coverage: use chr1 as a fast representative baseline
        ref_flag = ['--reference', cram_reference] if cram_reference else []
        cov = subprocess.run(
            ['samtools', 'coverage'] + ref_flag + ['-r', 'chr1', alignment_path],
            capture_output=True, text=True, env=env
        )
        if cov.returncode == 0:
            for line in cov.stdout.splitlines():
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 7:
                    try:
                        mean_coverage = float(parts[6])
                    except ValueError:
                        pass
                    break
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
