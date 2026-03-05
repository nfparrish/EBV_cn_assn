process HERITABILITY_ANALYSIS {
  label 'assoc_compute'
  tag 'heritability'

  input:
    path cn_normalized_table
    path pedigree_manifest

  output:
    path 'heritability_heritability_results.json', emit: results
    path 'heritability_heritability_trios.tsv', emit: trios
    path 'heritability_heritability_plot_data.json', emit: plot_data

  script:
  def cn_col = params.association_cn_column ?: 'depth_ratio_cn_norm'
  def n_perm = params.heritability_n_permutations ?: 10000
  def n_boot = params.heritability_n_bootstrap ?: 1000
  def seed = params.heritability_seed ?: 42
  """
  python3 ${projectDir}/scripts/compute_heritability.py \
    ${cn_normalized_table} \
    ${pedigree_manifest} \
    heritability \
    --cn-column '${cn_col}' \
    --n-permutations ${n_perm} \
    --n-bootstrap ${n_boot} \
    --seed ${seed}
  """
}
