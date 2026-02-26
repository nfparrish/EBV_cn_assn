nextflow.enable.dsl = 2

include { EBV_REFERENCE_PREP } from '../modules/local/ebv_reference_prep/main'
include { EBV_VARIANT_CALL } from '../modules/local/ebv_variant_call/main'
include { EBV_VARIANT_MATRIX } from '../modules/local/ebv_variant_matrix/main'
include { HOST_VARIANT_INGEST } from '../modules/local/host_variant_ingest/main'
include { HOST_VARIANT_HARMONIZE } from '../modules/local/host_variant_harmonize/main'
include { VARIANT_FILTER } from '../modules/local/variant_filter/main'
include { VARIANT_PROVENANCE } from '../modules/local/variant_provenance/main'

workflow EBV_VARIANT_LAYER {
  take:
    manifest_path

  main:
    if (!params.host_variant_input) {
      error "Missing required parameter for variant workflow: --host_variant_input <path/to/host_variant.tsv|vcf|plink-export>"
    }

    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    host_input_file = Channel.fromPath(params.host_variant_input, checkIfExists: true)
    manifest_rows = manifest_file.splitCsv(header: true)

    reference = EBV_REFERENCE_PREP(manifest_path)
    sample_variants = EBV_VARIANT_CALL(manifest_rows, reference.reference)
    matrix = EBV_VARIANT_MATRIX(sample_variants.variant_tsv.map { it[1] }.collect())
    host_ingest = HOST_VARIANT_INGEST(manifest_file, host_input_file)
    host_harmonized = HOST_VARIANT_HARMONIZE(manifest_file, host_ingest.staged_variants)
    filtered = VARIANT_FILTER(sample_variants.variant_tsv.map { it[1] }.collect(), host_harmonized.harmonized)
    provenance = VARIANT_PROVENANCE(
      filtered.counts,
      filtered.rules,
      matrix.summary,
      host_ingest.ingest_report,
      host_harmonized.harmonize_report,
    )

  emit:
    reference_manifest = reference.reference
    variant_catalog = matrix.variant_catalog
    variant_af_matrix = matrix.af_matrix
    variant_binary_matrix = matrix.binary_matrix
    variant_summary = matrix.summary
    host_variant_staged = host_ingest.staged_variants
    host_variant_sample_map = host_ingest.sample_map
    host_variant_ingest_report = host_ingest.ingest_report
    host_variant_harmonized = host_harmonized.harmonized
    host_variant_harmonize_report = host_harmonized.harmonize_report
    host_variant_sample_mismatch = host_harmonized.sample_mismatch
    ebv_variant_unfiltered = filtered.ebv_unfiltered
    ebv_variant_filtered = filtered.ebv_filtered
    host_variant_filtered = filtered.host_filtered
    variant_filter_rules = filtered.rules
    variant_filter_counts = filtered.counts
    variant_provenance_summary = provenance.summary
    variant_provenance_counts = provenance.counts
    variant_provenance_report = provenance.report
}
