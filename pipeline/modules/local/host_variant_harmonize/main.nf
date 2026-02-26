process HOST_VARIANT_HARMONIZE {
  label 'variant_compute'
  tag 'host_variant_harmonize'

  input:
    path manifest_file
    path staged_variants

  output:
    path 'host_variant_harmonized.tsv', emit: harmonized
    path 'host_variant_harmonize_report.json', emit: harmonize_report
    path 'host_variant_sample_mismatch.tsv', emit: sample_mismatch

  script:
  """
  python3 ${projectDir}/scripts/harmonize_host_variants.py ${manifest_file} ${staged_variants}
  """
}
