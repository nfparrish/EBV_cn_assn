process RNA_ASSOC_CN {
  label 'rna_compute'
  tag { "rna_assoc_cn_${normalization_label}" }

  input:
    path cn_input
    path rna_expression
    path sample_map
    val normalization_label

  output:
    path 'rna_assoc_cn_results.tsv', emit: results
    path 'rna_assoc_cn_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/run_rna_assoc_cn.py \
    ${cn_input} \
    ${rna_expression} \
    ${sample_map} \
    ${normalization_label} \
    rna_assoc_cn_results.tsv \
    rna_assoc_cn_report.json \
    ${params.rna_assoc_min_samples}
  """
}
