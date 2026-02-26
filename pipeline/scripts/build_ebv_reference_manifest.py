import json
import sys
from pathlib import Path

reference_accession = sys.argv[1]
viral_reference_fasta = sys.argv[2].strip() if len(sys.argv) > 2 else ''

if viral_reference_fasta:
    source_type = 'user_fasta'
    source = viral_reference_fasta
else:
    source_type = 'default_accession'
    source = reference_accession

payload = {
    'reference_accession': reference_accession,
    'source_type': source_type,
    'source': source,
}

Path('ebv_reference_manifest.json').write_text(json.dumps(payload, indent=2))
