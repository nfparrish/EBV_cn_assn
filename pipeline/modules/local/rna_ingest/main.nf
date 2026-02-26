process RNA_INGEST {
  label 'rna_compute'
  tag 'rna_ingest'

  input:
    path manifest_file

  output:
    path 'rna_expression_canonical.tsv', emit: expression_canonical
    path 'rna_sample_map.tsv', emit: sample_map
    path 'rna_ingest_report.json', emit: ingest_report

  script:
  """
  python3 ${projectDir}/scripts/build_rna_ingest.py \
    ${manifest_file} \
    '${params.rna_expression_input ?: ""}' \
    '${params.rna_metadata_input ?: ""}' \
    ${projectDir}
  """
}
