process HOST_VARIANT_INGEST {
  label 'variant_compute'
  tag 'host_variant_ingest'

  input:
    path manifest_file
    path host_input_file

  output:
    path 'host_variant_staged.tsv', emit: staged_variants
    path 'host_variant_sample_map.tsv', emit: sample_map
    path 'host_variant_ingest_report.json', emit: ingest_report

  script:
  """
  python3 ${projectDir}/scripts/stage_host_variant_input.py ${manifest_file} ${host_input_file} ${params.host_variant_format}
  """
}
