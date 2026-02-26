process ASSOC_RUN_HOST {
  label 'assoc_compute'
  tag 'assoc_run_host'

  input:
    path association_input
    path host_variant_filtered

  output:
    path 'assoc_host_results.tsv', emit: results
    path 'assoc_host_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/run_assoc_variants.py \
    ${association_input} \
    ${host_variant_filtered} \
    host \
    assoc_host_results.tsv \
    assoc_host_report.json \
    ${params.association_min_samples}
  """
}
