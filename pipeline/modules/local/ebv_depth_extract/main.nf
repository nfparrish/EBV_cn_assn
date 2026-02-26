process EBV_DEPTH_EXTRACT {
  label 'cn_compute'
  tag { sample_row.sample_id }

  input:
    val sample_row
    path reference_manifest

  output:
    tuple val(sample_row.sample_id), path("${sample_row.sample_id}.ebv_depth.json"), emit: depth

  script:
  """
  python3 ${projectDir}/scripts/compute_ebv_depth.py '${groovy.json.JsonOutput.toJson(sample_row)}' ${reference_manifest}
  """
}
