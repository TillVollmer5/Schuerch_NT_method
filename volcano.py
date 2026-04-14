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


# --- Class annotation helper --------------------------------------------------

def _load_class_annotation(cfg, feature_ids):
    """
    Build highlight_map and class_label_map from feature_metadata_enriched.csv.
    See pca.py for full documentation — this is an identical copy used here
    to avoid a cross-script import dependency.
    """
    highlight_rules = getattr(cfg, "CLASS_HIGHLIGHT", [])
    label_col       = getattr(cfg, "CLASS_LABEL_COLUMN", "")
    enriched_path   = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")

    highlight_map   = {}
    class_label_map = {}

    if not highlight_rules and not label_col:
        return highlight_map, class_label_map
    if not os.path.exists(enriched_path):
        return highlight_map, class_label_map

    enriched = pd.read_csv(enriched_path, index_col="feature_id")

    for rule in highlight_rules:
        col   = rule.get("column", "")
        val   = rule.get("value",  "")
        color = rule.get("color",  "#e74c3c")
        if col not in enriched.columns:
            continue
        for fid in feature_ids:
            if fid in enriched.index and enriched.loc[fid, col] == val:
                highlight_map[fid] = color

    if label_col and label_col in enriched.columns:
        for fid in feature_ids:
            if fid in enriched.index:
                v = enriched.loc[fid, label_col]
                if pd.notna(v) and str(v).strip():
                    class_label_map[fid] = str(v).strip()

    return highlight_map, class_label_map


# --- Feature label helper -----------------------------------------------------

def _load_feature_labels(cfg):
    """
    Load a feature_id -> display label mapping when FEATURE_LABEL = "name".

    Returns an empty dict when FEATURE_LABEL = "id" (the default) or when
    feature_name_map.csv is not found.  Individual entries fall back to the
    feature_id when no compound name is available.
    """
    if getattr(cfg, "FEATURE_LABEL", "id") != "name":
        return {}
    map_path = os.path.join(cfg.OUTPUT_DIR, "feature_name_map.csv")
    if not os.path.exists(map_path):
        print("  [info] feature_name_map.csv not found; using feature_id labels")
        return {}
    df = pd.read_csv(map_path, index_col="feature_id")
    if "compound_name" not in df.columns:
        return {}
    return {fid: name for fid, name in df["compound_name"].items()
            if pd.notna(name) and str(name).strip()}


# (normalization is now handled by normalization.py; volcano reads the
#  pre-built peak_matrix_processed_volcano.csv directly)


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


def _pvalue_for_feature(vals_a, vals_b, test):
    """
    Compute a raw p-value for a single feature given two sample arrays.

    test : str  — value of STAT_TEST_VOLCANO from config
      "mannwhitney"  - Mann-Whitney U (non-parametric; recommended default)
      "ttest"        - Welch's t-test (unequal variances)
      "ttest_equal"  - Student's t-test (equal variances assumed)
      "kruskal"      - Kruskal-Wallis (same as mannwhitney for two groups)
    """
    try:
        if test == "mannwhitney":
            _, p = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        elif test == "ttest":
            _, p = stats.ttest_ind(vals_a, vals_b, equal_var=False)
        elif test == "ttest_equal":
            _, p = stats.ttest_ind(vals_a, vals_b, equal_var=True)
        elif test == "kruskal":
            _, p = stats.kruskal(vals_a, vals_b)
        else:
            _, p = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
    except ValueError:
        p = 1.0
    return float(p)


def compute_volcano_stats(matrix_log2, samples_a, samples_b, stat_test="mannwhitney"):
    """
    Compute log2 fold change and per-feature p-values for a pairwise comparison.

    Parameters
    ----------
    matrix_log2 : DataFrame  features x samples  (log-normalized, no scaling)
    samples_a   : list of str  sample names for group A
    samples_b   : list of str  sample names for group B
    stat_test   : str  statistical test (see STAT_TEST_VOLCANO in config.py)

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
        pval   = _pvalue_for_feature(vals_a, vals_b, stat_test)

        records.append({"feature_id": fid, "log2FC": log2fc, "pvalue": pval})

    df               = pd.DataFrame(records).set_index("feature_id")
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
                 top_n, output_path, label_map=None,
                 highlight_map=None, class_label_map=None):
    """
    Draw and save a volcano plot.

    Parameters
    ----------
    results         : DataFrame  from classify(), indexed by feature_id
    group_a         : str        name of group A (positive log2FC direction)
    group_b         : str        name of group B
    fc_thresh       : float      log2 fold-change threshold
    p_thresh        : float      adjusted p-value threshold
    top_n           : int        number of top significant features to label
    output_path     : str
    label_map       : dict or None  feature_id -> display label (compound name)
    highlight_map   : dict or None  feature_id -> hex color; significant dots
                                    matching a rule get a colored ring outline
    class_label_map : dict or None  feature_id -> class string; appended as
                                    "[class]" to significant feature labels
    """
    highlight_map   = highlight_map   or {}
    class_label_map = class_label_map or {}

    x  = results["log2FC"].values
    # clip -log10(0) to a finite ceiling derived from the data
    adj_p      = results["adj_pvalue"].values.clip(1e-300, 1.0)
    y          = -np.log10(adj_p)
    directions = results["direction"].values
    feature_ids = results.index.tolist()
    # display labels: compound name when available, otherwise feature_id
    labels = ([label_map.get(fid, fid) for fid in feature_ids]
              if label_map else list(feature_ids))

    # append class suffix to labels where available
    if class_label_map:
        labels = [
            f"{lbl} [{class_label_map[fid]}]" if fid in class_label_map else lbl
            for lbl, fid in zip(labels, feature_ids)
        ]

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

    # colored ring outline on highlighted significant dots to show class
    if highlight_map:
        for i, fid in enumerate(feature_ids):
            if fid in highlight_map and directions[i] != "n.s.":
                ax.scatter(x[i], y[i], s=60, zorder=5,
                           facecolors="none",
                           edgecolors=highlight_map[fid],
                           linewidths=1.5)

    # threshold lines
    y_line = -np.log10(p_thresh)
    ax.axhline(y_line,    color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)
    ax.axvline( fc_thresh, color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)
    ax.axvline(-fc_thresh, color="#888888", linewidth=0.9,
               linestyle="--", zorder=2)

    # labels for top N significant features – rotated 270°, hanging below point
    sig_mask = directions != "n.s."
    if sig_mask.any() and top_n > 0:
        sig_idx    = np.where(sig_mask)[0]
        sig_sorted = sig_idx[np.argsort(adj_p[sig_idx])]
        label_idx  = sig_sorted[:top_n]

        ylim   = ax.get_ylim()
        y_drop = (ylim[1] - ylim[0]) * 0.04   # vertical offset below point
        xlim   = ax.get_xlim()
        x_step = (xlim[1] - xlim[0]) * 0.02  # horizontal nudge per collision

        # sort selected features left-to-right by x position so label order matches
        label_idx = label_idx[np.argsort(x[label_idx])]

        # initial x placements: same as point x
        placed = []   # list of (x_label, point_index)
        for i in label_idx:
            tx = x[i]
            # nudge right until no collision with already-placed labels
            for px, _ in placed:
                if abs(tx - px) < x_step:
                    tx = px + x_step
            placed.append((tx, i))

        for tx, i in placed:
            ty = y[i] - y_drop
            # text hangs below ty (va="top" with rotation 270 means text goes down)
            ax.text(tx, ty, labels[i],
                    fontsize=6, color="#333333", zorder=5,
                    rotation=270, ha="center", va="top")
            # connector line from point to top of label
            ax.plot([x[i], tx], [y[i], ty],
                    color="#aaaaaa", lw=0.5, zorder=4)

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
              loc="lower right")

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

    # Volcano reads the pre-built normalized matrix from normalization.py.
    # Normalization + log transform were applied there; scaling is always "none"
    # for volcano to preserve fold-change magnitudes.
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_volcano.csv")
    groups_path = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    for p in (matrix_path, groups_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} not found - run normalization.py first."
            )

    # samples x features (output of normalization.py)
    matrix_sf    = pd.read_csv(matrix_path, index_col="sample")
    group_series = pd.read_csv(groups_path,  index_col="sample")["group"]

    # align to samples present in both
    shared_samples = [s for s in matrix_sf.index if s in group_series.index]
    matrix_sf      = matrix_sf.loc[shared_samples]
    group_series   = group_series.reindex(shared_samples)

    # transpose to features x samples for compute_volcano_stats
    matrix_log2 = matrix_sf.T

    available_groups = sorted(group_series.unique())
    print(f"  features  : {matrix_log2.shape[0]}")
    print(f"  samples   : {matrix_log2.shape[1]}")
    print(f"  groups    : {available_groups}")

    # resolve comparisons
    comparisons = _resolve_comparisons(cfg.VOLCANO_COMPARISONS, available_groups)
    if not comparisons:
        print("  [warning] no valid comparisons to run.")
        return {}

    stat_test = getattr(cfg, "STAT_TEST_VOLCANO", "mannwhitney")
    print(f"  stat test  : {stat_test}  (BH-FDR corrected)")

    label_map   = _load_feature_labels(cfg) or None

    # class annotation (highlight + class label suffixes) — built once for all comparisons
    highlight_map, class_label_map = _load_class_annotation(cfg, matrix_log2.index)
    if highlight_map:
        print(f"  class highlight    : {len(highlight_map)} feature(s) highlighted")
    if class_label_map:
        print(f"  class labels       : {len(class_label_map)} feature(s) have class suffix")

    results_all = {}

    for group_a, group_b in comparisons:
        tag      = f"{group_a}_vs_{group_b}"
        samples_a = group_series[group_series == group_a].index.tolist()
        samples_b = group_series[group_series == group_b].index.tolist()

        print(f"\n  {group_a} (n={len(samples_a)}) vs {group_b} (n={len(samples_b)})")

        results = compute_volcano_stats(matrix_log2, samples_a, samples_b,
                                        stat_test=stat_test)
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
            group_a         = group_a,
            group_b         = group_b,
            fc_thresh       = cfg.VOLCANO_FC_THRESHOLD,
            p_thresh        = cfg.VOLCANO_P_THRESHOLD,
            top_n           = cfg.VOLCANO_TOP_LABELS,
            output_path     = out_png,
            label_map       = label_map,
            highlight_map   = highlight_map   or None,
            class_label_map = class_label_map or None,
        )
        print(f"    -> {out_png}")

        results_all[tag] = results

    return results_all


if __name__ == "__main__":
    run()
