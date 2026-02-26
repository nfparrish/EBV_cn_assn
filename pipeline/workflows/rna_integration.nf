nextflow.enable.dsl = 2

include { RNA_INGEST } from '../modules/local/rna_ingest/main'
include { RNA_COUNT_FASTQ } from '../modules/local/rna_count_fastq/main'
include { RNA_NORMALIZE } from '../modules/local/rna_normalize/main'
include { EBV_CN_ENGINE } from '../workflows/ebv_cn_engine'
include { EBV_REFERENCE_PREP } from '../modules/local/ebv_reference_prep/main'
include { EBV_VARIANT_CALL } from '../modules/local/ebv_variant_call/main'
include { RNA_ASSOC_CN as RNA_ASSOC_CN_LOG2 } from '../modules/local/rna_assoc_cn/main'
include { RNA_ASSOC_CN as RNA_ASSOC_CN_QUANT } from '../modules/local/rna_assoc_cn/main'
include { RNA_ASSOC_VARIANT as RNA_ASSOC_VARIANT_LOG2 } from '../modules/local/rna_assoc_variant/main'
include { RNA_ASSOC_VARIANT as RNA_ASSOC_VARIANT_QUANT } from '../modules/local/rna_assoc_variant/main'
include { RNA_ASSOC_QC } from '../modules/local/rna_assoc_qc/main'
include { RNA_VARIANT_PRESENCE_QC } from '../modules/local/rna_variant_presence_qc/main'

workflow RNA_INTEGRATION {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    manifest_rows = manifest_file.splitCsv(header: true)
    if (!(params.rna_input_mode in ['auto', 'matrix', 'fastq'])) {
      error "Unsupported --rna_input_mode '${params.rna_input_mode}'. Use 'auto', 'matrix', or 'fastq'."
    }
    use_fastq_counts = (params.rna_input_mode == 'fastq') || ((params.rna_input_mode == 'auto') && !params.rna_expression_input && params.rna_reference_fasta && params.rna_reference_gtf)
    cn_outputs = EBV_CN_ENGINE(manifest_path)
    reference = EBV_REFERENCE_PREP(manifest_path)
    sample_variants = EBV_VARIANT_CALL(manifest_rows, reference.reference)
    ebv_variant_inputs = sample_variants.variant_tsv.map { it[1] }.collect()
    expression_canonical = null
    sample_map = null
    ingest_report = null
    fastq_count_report = Channel.empty()

    if (use_fastq_counts) {
      ingest_fastq = RNA_COUNT_FASTQ(manifest_file)
      expression_canonical = ingest_fastq.expression_canonical
      sample_map = ingest_fastq.sample_map
      ingest_report = ingest_fastq.ingest_report
      fastq_count_report = ingest_fastq.count_report
    }
    else {
      ingest_matrix = RNA_INGEST(manifest_file)
      expression_canonical = ingest_matrix.expression_canonical
      sample_map = ingest_matrix.sample_map
      ingest_report = ingest_matrix.ingest_report
    }

    variant_presence = RNA_VARIANT_PRESENCE_QC(ebv_variant_inputs, sample_map)
    normalized = RNA_NORMALIZE(expression_canonical)
    cn_assoc_log2 = RNA_ASSOC_CN_LOG2(cn_outputs.cn_normalized_table, normalized.expression_log2, sample_map, 'log2_tpm_plus_1')
    cn_assoc_quant = RNA_ASSOC_CN_QUANT(cn_outputs.cn_normalized_table, normalized.expression_quantile, sample_map, 'quantile_normalized')
    variant_assoc_log2 = RNA_ASSOC_VARIANT_LOG2(ebv_variant_inputs, normalized.expression_log2, sample_map, 'log2_tpm_plus_1')
    variant_assoc_quant = RNA_ASSOC_VARIANT_QUANT(ebv_variant_inputs, normalized.expression_quantile, sample_map, 'quantile_normalized')
    qc = RNA_ASSOC_QC(
      cn_assoc_log2.results,
      cn_assoc_log2.report,
      cn_assoc_quant.results,
      cn_assoc_quant.report,
      variant_assoc_log2.results,
      variant_assoc_log2.report,
      variant_assoc_quant.results,
      variant_assoc_quant.report,
    )

  emit:
    rna_expression_canonical = expression_canonical
    rna_sample_map = sample_map
    rna_ingest_report = ingest_report
    rna_fastq_count_report = fastq_count_report
    rna_variant_presence_checks = variant_presence.checks
    rna_variant_presence_summary = variant_presence.summary
    rna_variant_presence_report = variant_presence.report
    rna_expression_log2tpm = normalized.expression_log2
    rna_expression_quantile = normalized.expression_quantile
    rna_normalization_report = normalized.normalization_report
    rna_assoc_cn_log2_results = cn_assoc_log2.results
    rna_assoc_cn_log2_report = cn_assoc_log2.report
    rna_assoc_cn_quantile_results = cn_assoc_quant.results
    rna_assoc_cn_quantile_report = cn_assoc_quant.report
    rna_assoc_variant_log2_results = variant_assoc_log2.results
    rna_assoc_variant_log2_report = variant_assoc_log2.report
    rna_assoc_variant_quantile_results = variant_assoc_quant.results
    rna_assoc_variant_quantile_report = variant_assoc_quant.report
    rna_assoc_qc_summary = qc.summary
    rna_assoc_qc_report = qc.report
    rna_assoc_qc_top_hits = qc.top_hits
}
