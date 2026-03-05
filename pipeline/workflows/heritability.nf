nextflow.enable.dsl = 2

include { EBV_CN_ENGINE } from '../workflows/ebv_cn_engine'
include { HERITABILITY_ANALYSIS } from '../modules/local/heritability_analysis/main'

workflow HERITABILITY {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)

    // Run CN engine to get normalized CN values
    cn_outputs = EBV_CN_ENGINE(manifest_path)

    // Run heritability analysis using CN table and pedigree manifest
    heritability = HERITABILITY_ANALYSIS(
      cn_outputs.cn_normalized_table,
      manifest_file,
    )

  emit:
    cn_normalized_table = cn_outputs.cn_normalized_table
    cn_summary = cn_outputs.cn_summary
    cn_qc_summary = cn_outputs.cn_qc_summary
    cn_qc_report = cn_outputs.cn_qc_report
    heritability_results = heritability.results
    heritability_trios = heritability.trios
    heritability_plot_data = heritability.plot_data
}
