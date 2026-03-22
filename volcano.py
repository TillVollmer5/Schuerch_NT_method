"""
volcano.py - Step 6 of the GCMS processing pipeline.

Produces a volcano plot for each requested pairwise group comparison.

For each comparison (group A vs group B):
  - log2 fold change  = mean(log2_A) - mean(log2_B)
    (difference of log2 means = log2 ratio of geometric means)
  - p-value           = Mann-Whitney U test per feature (non-parametric,
                        appropriate for small sample sizes)
  - adjusted p-value  = Benjamini-Hochberg FDR correction across all features

Input data: peak_matrix_blank_corrected.csv
  Sum normalization and log2(x+1) transformation are applied here before
  computing fold change and statistics.  The Pareto-scaled matrix is NOT used
  because scaling changes relative magnitudes between features, which would
  distort fold changes.

Feature classification:
  up in A    log2FC >  VOLCANO_FC_THRESHOLD  AND  adj-p < VOLCANO_P_THRESHOLD
  up in B    log2FC < -VOLCANO_FC_THRESHOLD  AND  adj-p < VOLCANO_P_THRESHOLD
  n.s.       all other features

The VOLCANO_TOP_LABELS most significant features (lowest adj-p) that exceed
both thresholds are labelled in the plot.

Output per comparison (A vs B):
  output/plots/volcano_A_vs_B.png   - volcano plot
  output/volcano_A_vs_B.csv         - full results table with columns:
      feature_id, log2FC, pvalue, adj_pvalue, direction

Input  : output/peak_matrix_blank_corrected.csv
         output/sample_groups.csv
Output : see above

Usage:
    python volcano.py
"""

import os
import sys
from itertools import combinations

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import config


# --- Normalization (applied inline, without Pareto scaling) ------------------

def _normalize_log2(matrix, normalization):
    """
    Apply sum or median normalization then log2(x+1) to a features x samples
    matrix.  Pareto scaling is intentionally omitted so that fold changes
    reflect real abundance ratios.
    """
    matrix = matrix.copy().apply(pd.to_numeric, errors="coerce").fillna(0)

    if normalization == "sum":
        col_sums   = matrix.sum(axis=0)
        median_sum = col_sums.median()
        if median_sum > 0:
            matrix = matrix.divide(col_sums.replace(0, np.nan), axis=1) * median_sum
            matrix = matrix.fillna(0)
    elif normalization == "median":
        col_medians   = matrix.median(axis=0)
        global_median = col_medians.median()
        if global_median > 0:
            col_medians = col_medians.replace(0, global_median)
            matrix      = matrix.divide(col_medians, axis=1) * global_median

    return np.log2(matrix + 1)


# --- Statistics ---------------------------------------------------------------

def _bh_correction(pvalues):
    """
    Benjamini-Hochberg FDR correction.
    Returns adjusted p-values clipped to [0, 1].
    """
    n = len(pvalues)
    if n == 0:
        return np.array([])
    order      = np.argsort(pvalues)
    ranks      = np.empty(n)
    ranks[order] = np.arange(1, n + 1)
    adj        = np.minimum(pvalues * n / ranks, 1.0)

    # enforce monotonicity across the sorted sequence
    adj_sorted = adj[order]
    for i in range(n - 2, -1, -1):
        if adj_sorted[i] > adj_sorted[i + 1]:
            adj_sorted[i] = adj_sorted[i + 1]
    adj[order] = adj_sorted
    return adj


def compute_volcano_stats(matrix_log2, samples_a, samples_b):
    """
    Compute log2 fold change and Mann-Whitney p-values for every feature.

    Parameters
    ----------
    matrix_log2 : DataFrame  features x samples  (log2-normalized)
    samples_a   : list of str  sample names for group A
    samples_b   : list of str  sample names for group B

    Returns
    -------
    DataFrame with columns: feature_id, log2FC, pvalue, adj_pvalue
    """
    a = matrix_log2[samples_a]
    b = matrix_log2[samples_b]

    records = []
    for fid in matrix_log2.index:
        vals_a = a.loc[fid].values.astype(float)
        vals_b = b.loc[fid].values.astype(float)

        log2fc = float(vals_a.mean() - vals_b.mean())

        # Mann-Whitney U; skip if either group has no variance (all zeros)
        try:
            _, pval = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        except ValueError:
            pval = 1.0

        records.append({"feature_id": fid, "log2FC": log2fc, "pvalue": float(pval)})

    df             = pd.DataFrame(records).set_index("feature_id")
    df["adj_pvalue"] = _bh_correction(df["pvalue"].values)
    return df


# --- Classification -----------------------------------------------------------

def classify(results, fc_thresh, p_thresh):
    """
    Add a 'direction' column: 'up_A', 'up_B', or 'n.s.'
    """
    up_a = (results["log2FC"]   >  fc_thresh) & (results["adj_pvalue"] < p_thresh)
    up_b = (results["log2FC"]   < -fc_thresh) & (results["adj_pvalue"] < p_thresh)

    results = results.copy()
    results["direction"] = "n.s."
    results.loc[up_a, "direction"] = "up_A"
    results.loc[up_b, "direction"] = "up_B"
    return results


# --- Plotting -----------------------------------------------------------------

_COL_UP_A = "#d6604d"   # red   - up in group A
_COL_UP_B = "#2166ac"   # blue  - up in group B
_COL_NS   = "#cccccc"   # grey  - not significant


def plot_volcano(results, group_a, group_b, fc_thresh, p_thresh,
                 top_n, output_path):
    """
    Draw and save a volcano plot.

    Parameters
    ----------
    results     : DataFrame  from classify(), indexed by feature_id
    group_a     : str        name of group A (positive log2FC direction)
    group_b     : str        name of group B
    fc_thresh   : float      log2 fold-change threshold
    p_thresh    : float      adjusted p-value threshold
    top_n       : int        number of top significant features to label
    output_path : str
    """
    x  = results["log2FC"].values
    # clip -log10(0) to a finite ceiling derived from the data
    adj_p      = results["adj_pvalue"].values.clip(1e-300, 1.0)
    y          = -np.log10(adj_p)
    directions = results["direction"].values
    feature_ids = results.index.tolist()

    colours = np.where(directions == "up_A", _COL_UP_A,
              np.where(directions == "up_B", _COL_UP_B, _COL_NS))
    sizes   = np.where(directions == "n.s.", 18, 28)
    zorders = np.where(directions == "n.s.", 1, 3)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    # plot n.s. first, then significant (so they render on top)
    for direction, col, sz, zo in [
        ("n.s.",  _COL_NS,   18, 1),
        ("up_A",  _COL_UP_A, 28, 3),
        ("up_B",  _COL_UP_B, 28, 3),
    ]:
        mask = directions == direction
        ax.scatter(x[mask], y[mask], c=col, s=sz, zorder=zo,
                   alpha=0.75, edgecolors="none")

    # threshold lines
    y_line = -np.log10(p_thresh)
    ax.axhline(y_line,    color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)
    ax.axvline( fc_thresh, color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)
    ax.axvline(-fc_thresh, color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)

    # labels for top N significant features
    sig_mask = directions != "n.s."
    if sig_mask.any() and top_n > 0:
        sig_idx   = np.where(sig_mask)[0]
        # sort significant features by adj_pvalue ascending (most significant first)
        sig_sorted = sig_idx[np.argsort(adj_p[sig_idx])]
        label_idx  = sig_sorted[:top_n]
        for i in label_idx:
            ax.annotate(
                feature_ids[i],
                (x[i], y[i]),
                textcoords="offset points",
                xytext=(5, 3),
                fontsize=6,
                color="#333333",
                zorder=5,
            )

    # counts in legend
    n_up_a = int((directions == "up_A").sum())
    n_up_b = int((directions == "up_B").sum())
    n_ns   = int((directions == "n.s.").sum())

    legend_handles = [
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=_COL_UP_A,
                      markersize=7,
                      label=f"up in {group_a}  (n={n_up_a})"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=_COL_UP_B,
                      markersize=7,
                      label=f"up in {group_b}  (n={n_up_b})"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor=_COL_NS,
                      markersize=7,
                      label=f"n.s.  (n={n_ns})"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, framealpha=0.9,
              loc="upper left")

    ax.set_xlabel(f"log2 fold change  ({group_a} / {group_b})", fontsize=10)
    ax.set_ylabel("-log10 (adj. p-value)", fontsize=10)
    ax.set_title(
        f"Volcano plot:  {group_a}  vs  {group_b}\n"
        f"FC threshold: {fc_thresh:.1f}   |   FDR threshold: {p_thresh}",
        fontsize=10, fontweight="bold",
    )
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.5)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# --- Comparison resolution ---------------------------------------------------

def _resolve_comparisons(comparisons, available_groups):
    """
    Resolve VOLCANO_COMPARISONS to a list of (groupA, groupB) tuples.

    Parameters
    ----------
    comparisons      : "all"  or  list of (groupA, groupB) tuples
    available_groups : list of str  groups present in the data
    """
    if comparisons == "all":
        return list(combinations(sorted(available_groups), 2))

    resolved = []
    for a, b in comparisons:
        missing = [g for g in (a, b) if g not in available_groups]
        if missing:
            print(f"  [skip] group(s) not found in data: {missing}  "
                  f"(requested comparison: {a} vs {b})")
            continue
        resolved.append((a, b))
    return resolved


# --- Main --------------------------------------------------------------------

def run(cfg=config):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 6: volcano plot ------------------------------------------")

    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    groups_path = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    for p in (matrix_path, groups_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} not found - run the preceding pipeline steps first."
            )

    # features x samples (raw blank-corrected areas)
    matrix       = pd.read_csv(matrix_path, index_col="feature_id")
    group_series = pd.read_csv(groups_path,  index_col="sample")["group"]

    # keep only samples present in the matrix
    shared_samples = [s for s in matrix.columns if s in group_series.index]
    matrix         = matrix[shared_samples]
    group_series   = group_series.reindex(shared_samples)

    available_groups = sorted(group_series.unique())
    print(f"  features  : {matrix.shape[0]}")
    print(f"  samples   : {matrix.shape[1]}")
    print(f"  groups    : {available_groups}")

    # apply normalization + log2 (no Pareto scaling)
    print(f"  normalization: {cfg.NORMALIZATION}  +  log2(x+1)")
    matrix_log2 = _normalize_log2(matrix, cfg.NORMALIZATION)

    # resolve comparisons
    comparisons = _resolve_comparisons(cfg.VOLCANO_COMPARISONS, available_groups)
    if not comparisons:
        print("  [warning] no valid comparisons to run.")
        return {}

    results_all = {}

    for group_a, group_b in comparisons:
        tag      = f"{group_a}_vs_{group_b}"
        samples_a = group_series[group_series == group_a].index.tolist()
        samples_b = group_series[group_series == group_b].index.tolist()

        print(f"\n  {group_a} (n={len(samples_a)}) vs {group_b} (n={len(samples_b)})")

        results = compute_volcano_stats(matrix_log2, samples_a, samples_b)
        results = classify(results,
                           fc_thresh=cfg.VOLCANO_FC_THRESHOLD,
                           p_thresh =cfg.VOLCANO_P_THRESHOLD)

        n_up_a = (results["direction"] == "up_A").sum()
        n_up_b = (results["direction"] == "up_B").sum()
        print(f"    up in {group_a}: {n_up_a}  |  up in {group_b}: {n_up_b}  "
              f"|  n.s.: {(results['direction'] == 'n.s.').sum()}")

        # save results table
        out_csv = os.path.join(cfg.OUTPUT_DIR, f"volcano_{tag}.csv")
        results.to_csv(out_csv)
        print(f"    -> {out_csv}")

        # save plot
        out_png = os.path.join(plots_dir, f"volcano_{tag}.png")
        plot_volcano(
            results,
            group_a      = group_a,
            group_b      = group_b,
            fc_thresh    = cfg.VOLCANO_FC_THRESHOLD,
            p_thresh     = cfg.VOLCANO_P_THRESHOLD,
            top_n        = cfg.VOLCANO_TOP_LABELS,
            output_path  = out_png,
        )
        print(f"    -> {out_png}")

        results_all[tag] = results

    return results_all


if __name__ == "__main__":
    run()
