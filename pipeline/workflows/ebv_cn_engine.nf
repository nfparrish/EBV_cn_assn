nextflow.enable.dsl = 2

include { EBV_REFERENCE_PREP } from '../modules/local/ebv_reference_prep/main'
include { EBV_DEPTH_EXTRACT } from '../modules/local/ebv_depth_extract/main'
include { EBV_CN_COMPUTE } from '../modules/local/ebv_cn_compute/main'
include { EBV_CN_NORMALIZE } from '../modules/local/ebv_cn_normalize/main'

workflow EBV_CN_ENGINE {
  take:
    manifest_path

  main:
    manifest_file = Channel.fromPath(manifest_path, checkIfExists: true)
    manifest_rows = manifest_file.splitCsv(header: true)

    reference = EBV_REFERENCE_PREP(manifest_path)
    depth = EBV_DEPTH_EXTRACT(manifest_rows, reference.reference)
    cn_raw = EBV_CN_COMPUTE(depth.depth)
    cn_norm = EBV_CN_NORMALIZE(cn_raw.cn_raw.map { it[1] }.collect())

  emit:
    reference_manifest = reference.reference
    cn_raw_table = cn_norm.raw_table
    cn_normalized_table = cn_norm.normalized_table
    cn_summary = cn_norm.summary
}
