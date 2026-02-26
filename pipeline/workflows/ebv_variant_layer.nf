nextflow.enable.dsl = 2

include { EBV_REFERENCE_PREP } from '../modules/local/ebv_reference_prep/main'
include { EBV_VARIANT_CALL } from '../modules/local/ebv_variant_call/main'
include { EBV_VARIANT_MATRIX } from '../modules/local/ebv_variant_matrix/main'

workflow EBV_VARIANT_LAYER {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    manifest_rows = manifest_file.splitCsv(header: true)

    reference = EBV_REFERENCE_PREP(manifest_path)
    sample_variants = EBV_VARIANT_CALL(manifest_rows, reference.reference)
    matrix = EBV_VARIANT_MATRIX(sample_variants.variant_tsv.map { it[1] }.collect())

  emit:
    reference_manifest = reference.reference
    variant_catalog = matrix.variant_catalog
    variant_af_matrix = matrix.af_matrix
    variant_binary_matrix = matrix.binary_matrix
    variant_summary = matrix.summary
}
