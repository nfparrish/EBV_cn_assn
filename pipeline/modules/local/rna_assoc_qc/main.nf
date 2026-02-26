process RNA_ASSOC_QC {
  label 'rna_compute'
  tag 'rna_assoc_qc'

  input:
    path cn_log2_results, stageAs: 'cn_log2_results.tsv'
    path cn_log2_report, stageAs: 'cn_log2_report.json'
    path cn_quant_results, stageAs: 'cn_quant_results.tsv'
    path cn_quant_report, stageAs: 'cn_quant_report.json'
    path var_log2_results, stageAs: 'var_log2_results.tsv'
    path var_log2_report, stageAs: 'var_log2_report.json'
    path var_quant_results, stageAs: 'var_quant_results.tsv'
    path var_quant_report, stageAs: 'var_quant_report.json'

  output:
    path 'rna_assoc_qc_summary.json', emit: summary
    path 'rna_assoc_qc_report.md', emit: report
    path 'rna_assoc_qc_top_hits.tsv', emit: top_hits

  script:
  """
  python3 ${projectDir}/scripts/build_rna_assoc_qc.py \
    ${cn_log2_results} ${cn_log2_report} \
    ${cn_quant_results} ${cn_quant_report} \
    ${var_log2_results} ${var_log2_report} \
    ${var_quant_results} ${var_quant_report}
  """
}
