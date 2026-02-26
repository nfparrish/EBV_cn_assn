process EBV_CN_COMPUTE {
  label 'cn_compute'
  tag { sample_id }

  input:
    tuple val(sample_id), path(depth_json)

  output:
    tuple val(sample_id), path("${sample_id}.ebv_cn_raw.json"), emit: cn_raw

  script:
  """
  python3 ${projectDir}/scripts/compute_ebv_cn_raw.py ${depth_json}
  """
}
