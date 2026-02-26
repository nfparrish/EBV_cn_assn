nextflow.enable.dsl = 2

include { EBV_CN_ENGINE } from '../workflows/ebv_cn_engine'
include { ASSOC_COVARIATE_HARMONIZE } from '../modules/local/assoc_covariate_harmonize/main'
include { ASSOC_CN_PREPARE } from '../modules/local/assoc_cn_prepare/main'

workflow WGS_ASSOCIATION {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    cn_outputs = EBV_CN_ENGINE(manifest_path)
    covariates = ASSOC_COVARIATE_HARMONIZE(manifest_file)
    prepared = ASSOC_CN_PREPARE(cn_outputs.cn_normalized_table, covariates.covariates)

  emit:
    cn_normalized_table = cn_outputs.cn_normalized_table
    cn_summary = cn_outputs.cn_summary
    association_covariates = covariates.covariates
    association_covariates_report = covariates.report
    association_input = prepared.model_input
    association_cn_prepare_report = prepared.report
}
