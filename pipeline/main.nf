nextflow.enable.dsl = 2

include { WGS_INTAKE } from './workflows/wgs_intake'
include { EBV_CN_ENGINE } from './workflows/ebv_cn_engine'
include { EBV_VARIANT_LAYER } from './workflows/ebv_variant_layer'
include { WGS_ASSOCIATION } from './workflows/wgs_association'

workflow {
  if (!params.manifest) {
    error "Missing required parameter: --manifest <path/to/manifest.csv>"
  }

  if (params.workflow == 'intake') {
    WGS_INTAKE(params.manifest)
  }
  else if (params.workflow == 'cn') {
    EBV_CN_ENGINE(params.manifest)
  }
  else if (params.workflow == 'variant') {
    EBV_VARIANT_LAYER(params.manifest)
  }
  else if (params.workflow == 'assoc') {
    WGS_ASSOCIATION(params.manifest)
  }
  else {
    error "Unsupported --workflow '${params.workflow}'. Use 'intake', 'cn', 'variant', or 'assoc'."
  }
}
