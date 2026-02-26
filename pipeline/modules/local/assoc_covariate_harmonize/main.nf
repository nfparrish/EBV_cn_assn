process ASSOC_COVARIATE_HARMONIZE {
  label 'assoc_compute'
  tag 'assoc_covariates'

  input:
    path manifest_file

  output:
    path 'association_covariates.tsv', emit: covariates
    path 'association_covariates_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/build_assoc_covariates.py \
    ${manifest_file} \
    '${params.association_covariates_input ?: ""}' \
    ${projectDir}
  """
}
