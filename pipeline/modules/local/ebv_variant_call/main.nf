process EBV_VARIANT_CALL {
  label 'variant_compute'
  tag { sample_row.sample_id }

  input:
    val sample_row
    path reference_manifest

  output:
    tuple val(sample_row.sample_id), path("${sample_row.sample_id}.ebv_variants.tsv"), emit: variant_tsv

  script:
  """
  python3 ${projectDir}/scripts/compute_ebv_variants_mutect2.py '${groovy.json.JsonOutput.toJson(sample_row)}' ${reference_manifest}
  """
}
