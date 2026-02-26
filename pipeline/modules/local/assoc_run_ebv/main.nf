process ASSOC_RUN_EBV {
  label 'assoc_compute'
  tag 'assoc_run_ebv'

  input:
    path association_input
    path ebv_variant_filtered

  output:
    path 'assoc_ebv_results.tsv', emit: results
    path 'assoc_ebv_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/run_assoc_variants.py \
    ${association_input} \
    ${ebv_variant_filtered} \
    ebv \
    assoc_ebv_results.tsv \
    assoc_ebv_report.json \
    ${params.association_min_samples}
  """
}
