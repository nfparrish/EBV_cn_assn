process VARIANT_FILTER {
  label 'variant_compute'
  tag 'variant_filter'

  input:
    path ebv_variant_tsvs
    path host_harmonized

  output:
    path 'ebv_variants_unfiltered.tsv', emit: ebv_unfiltered
    path 'ebv_variants_filtered.tsv', emit: ebv_filtered
    path 'host_variant_filtered.tsv', emit: host_filtered
    path 'variant_filter_rules.json', emit: rules
    path 'variant_filter_counts.json', emit: counts

  script:
  """
  python3 ${projectDir}/scripts/filter_variants.py \
    ${params.variant_filter_ebv_min_af} \
    ${params.variant_filter_ebv_min_dp} \
    ${params.variant_filter_ebv_min_qual} \
    ${params.variant_filter_host_min_dosage} \
    ${params.variant_filter_profile} \
    ${host_harmonized}
  """
}
