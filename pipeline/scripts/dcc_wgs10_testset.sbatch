#!/usr/bin/env bash
#SBATCH --job-name=gsd-wgs10
#SBATCH --partition=common
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR_DEFAULT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

REPO_DIR="${REPO_DIR:-${REPO_DIR_DEFAULT}}"
STUDY_ACCESSION="${STUDY_ACCESSION:-PRJNA41139}"
TESTSET_SIZE="${TESTSET_SIZE:-10}"

ACCOUNT="${ACCOUNT:-${SLURM_JOB_ACCOUNT:-}}"
PARTITION="${PARTITION:-${SLURM_JOB_PARTITION:-common}}"
QOS="${QOS:-}"

WORK_BASE="${WORK_BASE:-/work/${USER}/gsd_wgs10}"
DATA_DIR="${DATA_DIR:-${WORK_BASE}/data/wgs10}"
META_DIR="${META_DIR:-${WORK_BASE}/meta}"
NXF_WORK="${NXF_WORK:-${WORK_BASE}/nxf_work}"
SINGULARITY_CACHE="${SINGULARITY_CACHE:-${WORK_BASE}/singularity-cache}"
RESULTS_DIR="${RESULTS_DIR:-${WORK_BASE}/results/wgs10_report_$(date +%Y%m%d_%H%M%S)}"

MANIFEST_PATH="${MANIFEST_PATH:-${META_DIR}/wgs10_manifest.csv}"
HOST_VARIANTS_PATH="${HOST_VARIANTS_PATH:-${META_DIR}/wgs10_host_variants.tsv}"
DOWNLOAD_PLAN_PATH="${DOWNLOAD_PLAN_PATH:-${META_DIR}/wgs10_download_plan.tsv}"

if [[ -z "${ACCOUNT}" ]]; then
  echo "ERROR: ACCOUNT is not set. Submit with --account <your_account> or export ACCOUNT." >&2
  exit 1
fi

for cmd in nextflow java curl python3; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: required command not found: ${cmd}" >&2
    exit 1
  fi
done

mkdir -p "${DATA_DIR}" "${META_DIR}" "${NXF_WORK}" "${SINGULARITY_CACHE}" "${RESULTS_DIR}"

export NXF_HOME="${NXF_HOME:-${WORK_BASE}/.nextflow}"
mkdir -p "${NXF_HOME}"

echo "[INFO] Building ${TESTSET_SIZE}-sample WGS test-set plan from ENA study ${STUDY_ACCESSION}"

python3 - "${STUDY_ACCESSION}" "${TESTSET_SIZE}" "${DATA_DIR}" "${MANIFEST_PATH}" "${HOST_VARIANTS_PATH}" "${DOWNLOAD_PLAN_PATH}" <<'PY'
import csv
import hashlib
import io
import pathlib
import sys
import urllib.parse
import urllib.request

study = sys.argv[1]
target_n = int(sys.argv[2])
data_dir = pathlib.Path(sys.argv[3])
manifest_path = pathlib.Path(sys.argv[4])
host_variants_path = pathlib.Path(sys.argv[5])
download_plan_path = pathlib.Path(sys.argv[6])

query = (
    f'study_accession="{study}" AND '
    'library_source="GENOMIC" AND '
    'library_strategy="WGS" AND '
    'library_layout="PAIRED"'
)

params = {
    'result': 'read_run',
    'query': query,
    'fields': 'run_accession,sample_alias,fastq_ftp,fastq_bytes,study_accession',
    'format': 'tsv',
    'limit': '5000',
}

url = 'https://www.ebi.ac.uk/ena/portal/api/search?' + urllib.parse.urlencode(params)
with urllib.request.urlopen(url, timeout=120) as resp:
    body = resp.read().decode('utf-8', errors='replace')

reader = csv.DictReader(io.StringIO(body), delimiter='\t')
selected = []
seen_samples = set()

for row in reader:
    fastq_field = (row.get('fastq_ftp') or '').strip()
    if not fastq_field:
        continue

    files = [x.strip() for x in fastq_field.split(';') if x.strip()]
    read1 = next((x for x in files if x.endswith('_1.fastq.gz')), None)
    read2 = next((x for x in files if x.endswith('_2.fastq.gz')), None)
    if not read1 or not read2:
        continue

    sample_alias = (row.get('sample_alias') or '').strip()
    sample_id = sample_alias if sample_alias else row['run_accession']
    if sample_id in seen_samples:
        continue
    seen_samples.add(sample_id)

    run = row['run_accession']
    read1_url = f"https://{read1}" if not read1.startswith('http') else read1
    read2_url = f"https://{read2}" if not read2.startswith('http') else read2
    read1_path = data_dir / f"{run}_1.fastq.gz"
    read2_path = data_dir / f"{run}_2.fastq.gz"

    selected.append({
        'sample_id': sample_id,
        'dataset': 'OTHER_PUBLIC',
        'population': 'UNK',
        'sex': 'U',
        'sequence_lab': f'ENA_{study}',
        'read1_path': str(read1_path),
        'read2_path': str(read2_path),
        'bam_path': '',
        'md5': '',
        'run_accession': run,
        'sample_alias': sample_alias,
        'read1_url': read1_url,
        'read2_url': read2_url,
    })

    if len(selected) >= target_n:
        break

if len(selected) < target_n:
    raise SystemExit(
        f"Only found {len(selected)} paired WGS runs in study {study}; need {target_n}."
    )

manifest_cols = [
    'sample_id', 'dataset', 'population', 'sex', 'sequence_lab',
    'read1_path', 'read2_path', 'bam_path', 'md5',
]
with manifest_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=manifest_cols)
    writer.writeheader()
    for row in selected:
        writer.writerow({k: row[k] for k in manifest_cols})

plan_cols = ['sample_id', 'run_accession', 'sample_alias', 'read1_url', 'read2_url', 'read1_path', 'read2_path']
with download_plan_path.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=plan_cols, delimiter='\t')
    writer.writeheader()
    for row in selected:
        writer.writerow({k: row[k] for k in plan_cols})

with host_variants_path.open('w', newline='') as handle:
    handle.write('sample_id\tvariant_id\tchrom\tpos\tref\talt\tgt\tdosage\n')
    for row in selected:
        sample_id = row['sample_id']
        digest = hashlib.sha256(sample_id.encode('utf-8')).hexdigest()
        dosage1 = int(digest[0], 16) % 3
        dosage2 = int(digest[1], 16) % 3
        gt_map = {0: '0/0', 1: '0/1', 2: '1/1'}
        handle.write(f"{sample_id}\trsHost001\t1\t55516888\tA\tG\t{gt_map[dosage1]}\t{dosage1}\n")
        handle.write(f"{sample_id}\trsHost002\t2\t169991231\tC\tT\t{gt_map[dosage2]}\t{dosage2}\n")

print(f"Wrote manifest: {manifest_path}")
print(f"Wrote host variants: {host_variants_path}")
print(f"Wrote download plan: {download_plan_path}")
PY

echo "[INFO] Downloading FASTQs into ${DATA_DIR}"
tail -n +2 "${DOWNLOAD_PLAN_PATH}" | while IFS=$'\t' read -r sample_id run_accession sample_alias read1_url read2_url read1_path read2_path; do
  if [[ -s "${read1_path}" && -s "${read2_path}" ]]; then
    echo "[INFO] Reusing existing files for ${sample_id} (${run_accession})"
    continue
  fi

  echo "[INFO] Downloading ${sample_id} (${run_accession})"
  curl -L --fail --retry 4 --retry-delay 5 -o "${read1_path}" "${read1_url}"
  curl -L --fail --retry 4 --retry-delay 5 -o "${read2_path}" "${read2_url}"
done

echo "[INFO] Launching pipeline with dcc profile"
cd "${REPO_DIR}"

CMD=(
  nextflow run pipeline/main.nf
  --workflow report
  --manifest "${MANIFEST_PATH}"
  --host_variant_input "${HOST_VARIANTS_PATH}"
  --outdir "${RESULTS_DIR}"
  -profile dcc
  --slurm_account "${ACCOUNT}"
  --slurm_partition "${PARTITION}"
  --hpc_work_dir "${NXF_WORK}"
  --hpc_singularity_cache "${SINGULARITY_CACHE}"
  -with-report "${RESULTS_DIR}/nextflow-report.html"
  -with-timeline "${RESULTS_DIR}/nextflow-timeline.html"
  -with-trace "${RESULTS_DIR}/nextflow-trace.txt"
  -resume
)

if [[ -n "${QOS}" ]]; then
  CMD+=(--slurm_qos "${QOS}")
fi

printf '[INFO] Command: %q ' "${CMD[@]}"
printf '\n'
"${CMD[@]}"

echo "[INFO] Run complete. Results: ${RESULTS_DIR}"
