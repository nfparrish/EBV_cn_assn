import hashlib
import json
import pathlib
import sys
from datetime import datetime, timezone


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    files = [
        pathlib.Path(sys.argv[1]),
        pathlib.Path(sys.argv[2]),
        pathlib.Path(sys.argv[3]),
        pathlib.Path(sys.argv[4]),
        pathlib.Path(sys.argv[5]),
    ]

    rows = []
    for path in files:
        rows.append({
            'file_name': path.name,
            'size_bytes': path.stat().st_size,
            'sha256': sha256(path),
        })

    with pathlib.Path('report_reproducibility_checksums.tsv').open('w') as handle:
        handle.write('file_name\tsize_bytes\tsha256\n')
        for row in rows:
            handle.write(f"{row['file_name']}\t{row['size_bytes']}\t{row['sha256']}\n")

    manifest = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'artifact_count': len(rows),
        'artifacts': rows,
    }
    pathlib.Path('report_reproducibility_manifest.json').write_text(json.dumps(manifest, indent=2))

    report_lines = [
        '# Report Reproducibility Check',
        '',
        f"- artifacts_checked: {len(rows)}",
        '- checksum_algorithm: sha256',
        '- manifest_file: report_reproducibility_manifest.json',
        '- checksum_table: report_reproducibility_checksums.tsv',
        '',
        '## Checked artifacts',
    ]
    for row in rows:
        report_lines.append(f"- {row['file_name']} ({row['size_bytes']} bytes)")

    pathlib.Path('report_reproducibility_report.md').write_text('\n'.join(report_lines) + '\n')


if __name__ == '__main__':
    main()
