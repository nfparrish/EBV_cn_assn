#!/usr/bin/env python3
"""Estimate heritability of EBV copy number using parent-offspring trios.

Methods implemented (following Houldcroft et al., PLOS ONE 2014):
  1. Mid-parent regression: h² = slope(offspring_CN ~ mid_parent_CN)
     after residualizing for population, sex, and sequencing batch.
  2. Intraclass correlation coefficient (ICC) among family members.
  3. Permutation test: shuffle offspring labels N times to get null distribution.
  4. Bootstrap SE: resample trios with replacement for confidence intervals.

Usage:
    python compute_heritability.py <cn_table.tsv> <pedigree_manifest.csv> <output_prefix>
        [--cn-column depth_ratio_cn_norm]
        [--n-permutations 10000]
        [--n-bootstrap 1000]
        [--seed 42]

Inputs:
    cn_table.tsv         - Cohort CN table (from EBV_CN_ENGINE normalize step)
    pedigree_manifest.csv - Manifest with family_id, father_id, mother_id columns

Outputs:
    <prefix>_heritability_results.json - Point estimates, SE, p-values
    <prefix>_heritability_trios.tsv    - Per-trio CN values for plotting
    <prefix>_heritability_plot_data.json - Data for downstream figure generation
"""
import argparse
import csv
import json
import math
import pathlib
import random
import statistics
import sys


def read_tsv(path):
    """Read a TSV file into a list of dicts."""
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def read_csv(path):
    """Read a CSV file into a list of dicts."""
    with open(path) as fh:
        first_line = fh.readline()
        delimiter = "\t" if "\t" in first_line else ","
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        return list(reader)


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def residualize(values, groups):
    """Remove group means from values (simple categorical residualization).
    Used to adjust for population and batch effects before regression.
    """
    from collections import defaultdict
    group_sums = defaultdict(float)
    group_counts = defaultdict(int)
    for v, g in zip(values, groups):
        if v is not None:
            group_sums[g] += v
            group_counts[g] += 1
    group_means = {g: group_sums[g] / group_counts[g] for g in group_counts if group_counts[g] > 0}
    return [v - group_means.get(g, 0.0) if v is not None else None for v, g in zip(values, groups)]


def ols_slope_intercept(x, y):
    """Ordinary least squares: y = a + b*x. Returns (slope, intercept, r_squared)."""
    n = len(x)
    if n < 3:
        return float("nan"), float("nan"), float("nan")
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    ss_xy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    ss_xx = sum((xi - mean_x) ** 2 for xi in x)
    ss_yy = sum((yi - mean_y) ** 2 for yi in y)
    if ss_xx < 1e-15:
        return float("nan"), float("nan"), float("nan")
    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x
    r_squared = (ss_xy ** 2) / (ss_xx * ss_yy) if ss_yy > 1e-15 else float("nan")
    return slope, intercept, r_squared


def slope_se(x, y, slope, intercept):
    """Standard error of the OLS slope estimate."""
    n = len(x)
    if n < 3:
        return float("nan")
    residuals = [yi - (intercept + slope * xi) for xi, yi in zip(x, y)]
    mse = sum(r ** 2 for r in residuals) / (n - 2)
    mean_x = sum(x) / n
    ss_xx = sum((xi - mean_x) ** 2 for xi in x)
    if ss_xx < 1e-15:
        return float("nan")
    return math.sqrt(mse / ss_xx)


def compute_icc(family_values):
    """Compute intraclass correlation coefficient (ICC(1,1)) using one-way ANOVA.
    family_values: dict of family_id -> list of CN values for family members.
    """
    k_groups = len(family_values)
    if k_groups < 2:
        return float("nan"), float("nan")

    all_values = []
    group_means = {}
    for fid, vals in family_values.items():
        if len(vals) < 2:
            continue
        group_means[fid] = sum(vals) / len(vals)
        all_values.extend(vals)

    if not all_values or len(group_means) < 2:
        return float("nan"), float("nan")

    grand_mean = sum(all_values) / len(all_values)
    n_total = len(all_values)
    k = len(group_means)

    # Between-group sum of squares
    ss_between = sum(
        len(family_values[fid]) * (gm - grand_mean) ** 2
        for fid, gm in group_means.items()
    )
    # Within-group sum of squares
    ss_within = sum(
        (v - group_means.get(fid, grand_mean)) ** 2
        for fid, vals in family_values.items()
        if fid in group_means
        for v in vals
    )

    df_between = k - 1
    df_within = n_total - k
    if df_between < 1 or df_within < 1:
        return float("nan"), float("nan")

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    # Average group size
    n_sizes = [len(family_values[fid]) for fid in group_means]
    n_0 = (n_total - sum(ni ** 2 for ni in n_sizes) / n_total) / (k - 1)

    if n_0 < 1e-10:
        return float("nan"), float("nan")

    icc = (ms_between - ms_within) / (ms_between + (n_0 - 1) * ms_within)

    # F-test p-value
    if ms_within > 0:
        f_stat = ms_between / ms_within
        p_value = f_test_pvalue(f_stat, df_between, df_within)
    else:
        p_value = float("nan")

    return icc, p_value


def f_test_pvalue(f_stat, df1, df2):
    """Approximate p-value for F-test using Wilson-Hilferty normal approximation."""
    if f_stat <= 0 or df1 <= 0 or df2 <= 0:
        return float("nan")
    nu1, nu2 = float(df1), float(df2)
    z = (
        ((f_stat / nu1 * nu2) ** (1.0 / 3) - (1 - 2.0 / (9 * nu2)))
        / math.sqrt(2.0 / (9 * nu2))
    )
    p = 0.5 * math.erfc(z / math.sqrt(2))
    return max(0.0, min(1.0, p))


def main():
    parser = argparse.ArgumentParser(description="EBV CN heritability estimation")
    parser.add_argument("cn_table", help="Cohort CN table (TSV)")
    parser.add_argument("pedigree_manifest", help="Manifest with pedigree columns (CSV)")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("--cn-column", default="depth_ratio_cn_norm")
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    # --- Load data ---
    cn_rows = read_tsv(args.cn_table)
    manifest_rows = read_csv(args.pedigree_manifest)

    cn_by_sample = {}
    for row in cn_rows:
        sid = (row.get("sample_id") or "").strip()
        val = to_float(row.get(args.cn_column))
        if sid and val is not None:
            cn_by_sample[sid] = {
                "cn": val,
                "population": row.get("population", ""),
                "sex": row.get("sex", ""),
                "sequence_lab": row.get("sequence_lab", ""),
            }

    # Build pedigree from manifest
    manifest_by_id = {}
    for row in manifest_rows:
        sid = (row.get("sample_id") or "").strip()
        if sid:
            manifest_by_id[sid] = row

    # --- Identify complete trios with CN data ---
    trios = []
    for sid, mrow in manifest_by_id.items():
        father_id = (mrow.get("father_id") or "").strip()
        mother_id = (mrow.get("mother_id") or "").strip()
        family_id = (mrow.get("family_id") or "").strip()

        if not father_id or father_id == "0" or not mother_id or mother_id == "0":
            continue
        if sid not in cn_by_sample or father_id not in cn_by_sample or mother_id not in cn_by_sample:
            continue

        trios.append({
            "family_id": family_id,
            "child_id": sid,
            "father_id": father_id,
            "mother_id": mother_id,
            "child_cn": cn_by_sample[sid]["cn"],
            "father_cn": cn_by_sample[father_id]["cn"],
            "mother_cn": cn_by_sample[mother_id]["cn"],
            "mid_parent_cn": (cn_by_sample[father_id]["cn"] + cn_by_sample[mother_id]["cn"]) / 2.0,
            "child_pop": cn_by_sample[sid]["population"],
            "child_sex": cn_by_sample[sid]["sex"],
            "child_lab": cn_by_sample[sid]["sequence_lab"],
        })

    n_trios = len(trios)
    print(f"Found {n_trios} complete trios with CN data")

    if n_trios < 3:
        results = {
            "status": "insufficient_data",
            "n_trios": n_trios,
            "message": "Need at least 3 complete trios for heritability estimation",
        }
        pathlib.Path(f"{args.output_prefix}_heritability_results.json").write_text(
            json.dumps(results, indent=2)
        )
        return

    # --- Residualize CN for population + batch ---
    all_sample_ids = set()
    for t in trios:
        all_sample_ids.update([t["child_id"], t["father_id"], t["mother_id"]])

    group_key = {}
    for sid in all_sample_ids:
        if sid in cn_by_sample:
            group_key[sid] = f"{cn_by_sample[sid]['population']}_{cn_by_sample[sid]['sequence_lab']}"

    sample_list = sorted(all_sample_ids)
    raw_cn = [cn_by_sample[sid]["cn"] for sid in sample_list]
    groups = [group_key.get(sid, "unknown") for sid in sample_list]
    resid_cn = residualize(raw_cn, groups)
    resid_by_sample = dict(zip(sample_list, resid_cn))

    # --- Mid-parent regression (Falconer's method) ---
    mid_parent = [(t["father_cn"] + t["mother_cn"]) / 2.0 for t in trios]
    offspring = [t["child_cn"] for t in trios]

    mid_parent_resid = [
        (resid_by_sample.get(t["father_id"], 0) + resid_by_sample.get(t["mother_id"], 0)) / 2.0
        for t in trios
    ]
    offspring_resid = [resid_by_sample.get(t["child_id"], 0) for t in trios]

    h2_raw, intercept_raw, r2_raw = ols_slope_intercept(mid_parent, offspring)
    se_raw = slope_se(mid_parent, offspring, h2_raw, intercept_raw)

    h2_adj, intercept_adj, r2_adj = ols_slope_intercept(mid_parent_resid, offspring_resid)
    se_adj = slope_se(mid_parent_resid, offspring_resid, h2_adj, intercept_adj)

    def slope_pvalue(slope, se, n):
        if se <= 0 or math.isnan(se) or n < 3:
            return float("nan")
        t_stat = slope / se
        return math.erfc(abs(t_stat) / math.sqrt(2))

    p_raw = slope_pvalue(h2_raw, se_raw, n_trios)
    p_adj = slope_pvalue(h2_adj, se_adj, n_trios)

    # --- Permutation test ---
    observed_slope = h2_adj
    perm_slopes = []
    offspring_resid_copy = list(offspring_resid)
    for _ in range(args.n_permutations):
        random.shuffle(offspring_resid_copy)
        s, _, _ = ols_slope_intercept(mid_parent_resid, offspring_resid_copy)
        perm_slopes.append(s)

    perm_p = sum(1 for s in perm_slopes if s >= observed_slope) / args.n_permutations

    # --- Bootstrap SE and CI ---
    boot_slopes = []
    indices = list(range(n_trios))
    for _ in range(args.n_bootstrap):
        sample = random.choices(indices, k=n_trios)
        x_boot = [mid_parent_resid[i] for i in sample]
        y_boot = [offspring_resid[i] for i in sample]
        s, _, _ = ols_slope_intercept(x_boot, y_boot)
        if not math.isnan(s):
            boot_slopes.append(s)

    boot_se = statistics.stdev(boot_slopes) if len(boot_slopes) > 1 else float("nan")
    boot_slopes_sorted = sorted(boot_slopes)
    ci_low = boot_slopes_sorted[int(0.025 * len(boot_slopes_sorted))] if boot_slopes_sorted else float("nan")
    ci_high = boot_slopes_sorted[int(0.975 * len(boot_slopes_sorted))] if boot_slopes_sorted else float("nan")

    # --- ICC ---
    family_cn_values = {}
    for t in trios:
        fid = t["family_id"]
        if fid not in family_cn_values:
            family_cn_values[fid] = []
        family_cn_values[fid].extend([
            resid_by_sample.get(t["child_id"], 0),
            resid_by_sample.get(t["father_id"], 0),
            resid_by_sample.get(t["mother_id"], 0),
        ])

    icc_value, icc_pvalue = compute_icc(family_cn_values)

    # --- Write results ---
    results = {
        "status": "success",
        "n_trios": n_trios,
        "n_unique_samples": len(all_sample_ids),
        "cn_column": args.cn_column,
        "midparent_regression_raw": {
            "h2_estimate": round(h2_raw, 6) if not math.isnan(h2_raw) else None,
            "se": round(se_raw, 6) if not math.isnan(se_raw) else None,
            "r_squared": round(r2_raw, 6) if not math.isnan(r2_raw) else None,
            "p_value": round(p_raw, 8) if not math.isnan(p_raw) else None,
            "intercept": round(intercept_raw, 6) if not math.isnan(intercept_raw) else None,
        },
        "midparent_regression_adjusted": {
            "h2_estimate": round(h2_adj, 6) if not math.isnan(h2_adj) else None,
            "se_ols": round(se_adj, 6) if not math.isnan(se_adj) else None,
            "se_bootstrap": round(boot_se, 6) if not math.isnan(boot_se) else None,
            "ci_95_low": round(ci_low, 6) if not math.isnan(ci_low) else None,
            "ci_95_high": round(ci_high, 6) if not math.isnan(ci_high) else None,
            "r_squared": round(r2_adj, 6) if not math.isnan(r2_adj) else None,
            "p_value_ols": round(p_adj, 8) if not math.isnan(p_adj) else None,
            "p_value_permutation": round(perm_p, 8),
            "n_permutations": args.n_permutations,
            "n_bootstrap": args.n_bootstrap,
            "intercept": round(intercept_adj, 6) if not math.isnan(intercept_adj) else None,
            "adjustment": "population_x_batch_residualization",
        },
        "icc": {
            "icc_1_1": round(icc_value, 6) if not math.isnan(icc_value) else None,
            "p_value": round(icc_pvalue, 8) if not math.isnan(icc_pvalue) else None,
            "n_families": len(family_cn_values),
        },
        "method_notes": {
            "falconer": "h² estimated as slope of offspring CN on mid-parent CN",
            "adjustment": "CN residualized for population × sequencing batch before regression",
            "reference": "Houldcroft et al. PLOS ONE 2014; Mandage et al. PLOS ONE 2017",
            "icc": "ICC(1,1) from one-way random-effects ANOVA across families",
            "permutation": "Offspring labels shuffled to generate null distribution of slopes",
            "bootstrap": "Trios resampled with replacement for SE and 95% CI",
        },
        "seed": args.seed,
    }

    out_results = pathlib.Path(f"{args.output_prefix}_heritability_results.json")
    out_results.write_text(json.dumps(results, indent=2))
    print(f"Results written to {out_results}")

    # --- Write per-trio table ---
    trio_header = [
        "family_id", "child_id", "father_id", "mother_id",
        "child_cn", "father_cn", "mother_cn", "mid_parent_cn",
        "child_cn_resid", "father_cn_resid", "mother_cn_resid", "mid_parent_cn_resid",
        "child_pop", "child_sex", "child_lab",
    ]
    out_trios = pathlib.Path(f"{args.output_prefix}_heritability_trios.tsv")
    with out_trios.open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=trio_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, t in enumerate(trios):
            writer.writerow({
                "family_id": t["family_id"],
                "child_id": t["child_id"],
                "father_id": t["father_id"],
                "mother_id": t["mother_id"],
                "child_cn": f"{t['child_cn']:.8f}",
                "father_cn": f"{t['father_cn']:.8f}",
                "mother_cn": f"{t['mother_cn']:.8f}",
                "mid_parent_cn": f"{t['mid_parent_cn']:.8f}",
                "child_cn_resid": f"{offspring_resid[i]:.8f}",
                "father_cn_resid": f"{resid_by_sample.get(t['father_id'], 0):.8f}",
                "mother_cn_resid": f"{resid_by_sample.get(t['mother_id'], 0):.8f}",
                "mid_parent_cn_resid": f"{mid_parent_resid[i]:.8f}",
                "child_pop": t["child_pop"],
                "child_sex": t["child_sex"],
                "child_lab": t["child_lab"],
            })

    # --- Write plot data ---
    plot_data = {
        "midparent_regression": {
            "x": [round(v, 6) for v in mid_parent_resid],
            "y": [round(v, 6) for v in offspring_resid],
            "slope": round(h2_adj, 6) if not math.isnan(h2_adj) else None,
            "intercept": round(intercept_adj, 6) if not math.isnan(intercept_adj) else None,
            "labels": [t["child_id"] for t in trios],
            "populations": [t["child_pop"] for t in trios],
        },
        "permutation_distribution": {
            "observed_slope": round(observed_slope, 6) if not math.isnan(observed_slope) else None,
            "null_slopes_summary": {
                "mean": round(statistics.mean(perm_slopes), 6) if perm_slopes else None,
                "sd": round(statistics.stdev(perm_slopes), 6) if len(perm_slopes) > 1 else None,
                "percentiles": {
                    "2.5": round(sorted(perm_slopes)[int(0.025 * len(perm_slopes))], 6) if perm_slopes else None,
                    "97.5": round(sorted(perm_slopes)[int(0.975 * len(perm_slopes))], 6) if perm_slopes else None,
                },
            },
            "n_permutations": args.n_permutations,
            "p_value": round(perm_p, 8),
        },
        "bootstrap_distribution": {
            "slopes_summary": {
                "mean": round(statistics.mean(boot_slopes), 6) if boot_slopes else None,
                "se": round(boot_se, 6) if not math.isnan(boot_se) else None,
                "ci_95": [round(ci_low, 6) if not math.isnan(ci_low) else None,
                          round(ci_high, 6) if not math.isnan(ci_high) else None],
            },
            "n_bootstrap": args.n_bootstrap,
        },
    }
    out_plot = pathlib.Path(f"{args.output_prefix}_heritability_plot_data.json")
    out_plot.write_text(json.dumps(plot_data, indent=2))
    print(f"Plot data written to {out_plot}")


if __name__ == "__main__":
    main()
