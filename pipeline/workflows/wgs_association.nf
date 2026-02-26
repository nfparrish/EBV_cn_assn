nextflow.enable.dsl = 2

include { EBV_CN_ENGINE } from '../workflows/ebv_cn_engine'
include { EBV_VARIANT_LAYER } from '../workflows/ebv_variant_layer'
include { ASSOC_COVARIATE_HARMONIZE } from '../modules/local/assoc_covariate_harmonize/main'
include { ASSOC_CN_PREPARE } from '../modules/local/assoc_cn_prepare/main'
include { ASSOC_RUN_EBV } from '../modules/local/assoc_run_ebv/main'
include { ASSOC_RUN_HOST } from '../modules/local/assoc_run_host/main'

workflow WGS_ASSOCIATION {
  take:
    manifest_path

  main:
    if (!params.host_variant_input) {
      error "Missing required parameter for association workflow: --host_variant_input <path/to/host_variant.tsv|vcf|plink-export>"
    }

    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    cn_outputs = EBV_CN_ENGINE(manifest_path)
    variant_outputs = EBV_VARIANT_LAYER(manifest_path)
    covariates = ASSOC_COVARIATE_HARMONIZE(manifest_file)
    prepared = ASSOC_CN_PREPARE(cn_outputs.cn_normalized_table, covariates.covariates)
    ebv_assoc = ASSOC_RUN_EBV(prepared.model_input, variant_outputs.ebv_variant_filtered)
    host_assoc = ASSOC_RUN_HOST(prepared.model_input, variant_outputs.host_variant_filtered)

  emit:
    cn_normalized_table = cn_outputs.cn_normalized_table
    cn_summary = cn_outputs.cn_summary
    association_covariates = covariates.covariates
    association_covariates_report = covariates.report
    association_input = prepared.model_input
    association_cn_prepare_report = prepared.report
    assoc_ebv_results = ebv_assoc.results
    assoc_ebv_report = ebv_assoc.report
    assoc_host_results = host_assoc.results
    assoc_host_report = host_assoc.report
}
