nextflow.enable.dsl = 2

include { RNA_INGEST } from '../modules/local/rna_ingest/main'
include { RNA_NORMALIZE } from '../modules/local/rna_normalize/main'

workflow RNA_INTEGRATION {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    ingest = RNA_INGEST(manifest_file)
    normalized = RNA_NORMALIZE(ingest.expression_canonical)

  emit:
    rna_expression_canonical = ingest.expression_canonical
    rna_sample_map = ingest.sample_map
    rna_ingest_report = ingest.ingest_report
    rna_expression_log2tpm = normalized.expression_log2
    rna_expression_quantile = normalized.expression_quantile
    rna_normalization_report = normalized.normalization_report
}
