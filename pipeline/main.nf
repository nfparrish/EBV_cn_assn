nextflow.enable.dsl = 2

include { WGS_INTAKE } from './workflows/wgs_intake'

workflow {
  if (!params.manifest) {
    error "Missing required parameter: --manifest <path/to/manifest.csv>"
  }

  WGS_INTAKE(params.manifest)
}
