process QC_SAMPLE {
  label 'staging'
  tag { sample_id }

  input:
    tuple val(sample_id), path(staged_json), path(status_json)

  output:
    tuple val(sample_id), path("${sample_id}.qc.json"), emit: qc

  script:
  def cram_ref = params.cram_reference_fasta ? params.cram_reference_fasta : ''
  """
  python3 ${projectDir}/scripts/build_qc_sample.py \
    '${sample_id}' \
    '${staged_json}' \
    '${status_json}' \
    '${cram_ref}'
  """
}
