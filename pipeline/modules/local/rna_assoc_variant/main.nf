process RNA_ASSOC_VARIANT {
  label 'rna_compute'
  tag { "rna_assoc_variant_${normalization_label}" }

  input:
    path variant_inputs
    path rna_expression
    path sample_map
    val normalization_label

  output:
    path 'rna_assoc_variant_results.tsv', emit: results
    path 'rna_assoc_variant_report.json', emit: report

  script:
  """
  python3 ${projectDir}/scripts/run_rna_assoc_variant.py \
    ${rna_expression} \
    ${sample_map} \
    ${normalization_label} \
    rna_assoc_variant_results.tsv \
    rna_assoc_variant_report.json \
    ${params.rna_assoc_min_samples} \
    ${params.variant_filter_ebv_min_af} \
    ${params.variant_filter_ebv_min_dp} \
    ${params.variant_filter_ebv_min_qual}
  """
}
