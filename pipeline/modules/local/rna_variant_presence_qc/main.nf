process RNA_VARIANT_PRESENCE_QC {
  label 'rna_compute'
  tag 'rna_variant_presence_qc'

  input:
    path sample_variant_tsvs
    path sample_map

  output:
    path 'rna_variant_presence_checks.tsv', emit: checks
    path 'rna_variant_presence_summary.json', emit: summary
    path 'rna_variant_presence_report.md', emit: report

  script:
  """
  python3 ${projectDir}/scripts/check_rna_variant_presence.py \
    ${sample_map} \
    '${params.rna_expected_coding_variants ?: ""}' \
    '${params.rna_expected_donor_variants_input ?: ""}' \
    ${projectDir} \
    ${params.variant_filter_ebv_min_af} \
    ${params.variant_filter_ebv_min_dp} \
    ${params.variant_filter_ebv_min_qual}
  """
}
