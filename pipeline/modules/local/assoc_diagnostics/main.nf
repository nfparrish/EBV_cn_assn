process ASSOC_DIAGNOSTICS {
  label 'assoc_compute'
  tag 'assoc_diagnostics'

  input:
    path ebv_results
    path host_results

  output:
    path 'association_diagnostics.tsv', emit: diagnostics_tsv
    path 'association_diagnostics_summary.json', emit: diagnostics_summary

  script:
  """
  python3 ${projectDir}/scripts/build_assoc_diagnostics.py ${ebv_results} ${host_results}
  """
}
