process EBV_REFERENCE_PREP {
  label 'cn_compute'
  tag 'ebv_reference'

  input:
    val manifest_path

  output:
    path 'ebv_reference_manifest.json', emit: reference

  script:
  """
  python3 ${projectDir}/scripts/build_ebv_reference_manifest.py '${params.reference_accession}' '${params.viral_reference_fasta ?: ""}'
  """
}
