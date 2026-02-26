process RNA_COUNT_FASTQ {
  label 'rna_compute'
  tag 'rna_count_fastq'

  input:
    path manifest_file

  output:
    path 'rna_expression_canonical.tsv', emit: expression_canonical
    path 'rna_sample_map.tsv', emit: sample_map
    path 'rna_ingest_report.json', emit: ingest_report
    path 'rna_fastq_count_report.json', emit: count_report

  script:
  """
  python3 ${projectDir}/scripts/build_rna_counts_from_fastq.py \
    ${manifest_file} \
    '${params.rna_metadata_input ?: ""}' \
    '${params.rna_reference_fasta ?: ""}' \
    '${params.rna_reference_gtf ?: ""}' \
    ${projectDir} \
    ${params.rna_require_ebv_genes} \
    ${params.rna_threads} \
    '${params.rna_aligner ?: "minimap2"}' \
    '${params.rna_counter ?: "featureCounts"}'
  """
}
