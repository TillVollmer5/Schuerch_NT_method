"""
hca.py - Step 5 of the GCMS processing pipeline.

Performs Hierarchical Cluster Analysis (HCA) on the processed peak matrix and
writes a bidirectional clustered heatmap with dendrograms on both sample and
feature axes:

  Rows    = samples, annotated by group colour
  Columns = features; co-varying features cluster together

Input data: peak_matrix_processed.csv
  Normalized, log2-transformed, and Pareto-scaled (Step 3 output).
  This scaling is required - do NOT use the raw or blank-corrected matrix.

Output:
  output/plots/hca_heatmap.png   - clustered heatmap with bidirectional dendrograms
  output/hca_row_order.csv       - sample order from the row dendrogram
  output/hca_col_order.csv       - feature order from the column dendrogram

The dendrogram order CSVs are useful for:
  - Selecting co-varying metabolite clusters for Spearman correlation analysis
  - Reporting the exact cluster order in supplementary tables

Algorithm: Ward linkage + Euclidean distance by default (configurable via
config.py). Note: Ward linkage requires Euclidean distance; use "average" or
"complete" linkage if you want to switch to correlation distance.

Input  : output/peak_matrix_processed.csv
         output/sample_groups.csv
Output : see above

Usage:
    python hca.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive backend - no display needed
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

import config


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


# --- Colour palette (matches pca.py) -----------------------------------------

_PALETTE = [
    "#2166ac",   # blue
    "#d6604d",   # red
    "#4dac26",   # green
    "#8073ac",   # purple
    "#f4a582",   # orange
    "#1a1a1a",   # black
]


def _group_colours(group_series):
    """
    Map each group label to a colour from _PALETTE.

    Returns
    -------
    colour_series : Series  sample -> hex colour  (for seaborn row_colors)
    colour_map    : dict    group  -> hex colour  (for legend)
    """
    unique     = list(dict.fromkeys(group_series))
    colour_map = {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(unique)}
    return group_series.map(colour_map), colour_map


# --- Validation ---------------------------------------------------------------

def _validate_linkage_metric(linkage, metric):
    """
    Warn when an incompatible linkage / metric combination is requested.
    Ward linkage only supports Euclidean distance.
    """
    if linkage == "ward" and metric != "euclidean":
        raise ValueError(
            f"Ward linkage requires metric='euclidean', got '{metric}'. "
            "Change HCA_LINKAGE to 'average' or 'complete' in config.py "
            "to use a different distance metric."
        )


# --- HCA / clustermap --------------------------------------------------------

def plot_hca(matrix, group_series, linkage, metric, cmap,
             max_feature_labels, output_path):
    """
    Draw and save a bidirectional clustered heatmap with group row annotation.

    Parameters
    ----------
    matrix             : DataFrame  samples x features  (processed, scaled)
    group_series       : Series     sample -> group label
    linkage            : str        scipy linkage method: "ward", "average",
                                    "complete", "single"
    metric             : str        distance metric: "euclidean", "correlation"
    cmap               : str        matplotlib/seaborn diverging colormap
    max_feature_labels : int        label column axis when n_features <= this;
                                    0 = never show feature labels
    output_path        : str

    Returns
    -------
    g : seaborn.matrix.ClusterGrid
        The grid object; use g.dendrogram_row / g.dendrogram_col to extract
        reordered indices.
    """
    _validate_linkage_metric(linkage, metric)

    row_colours, colour_map = _group_colours(
        group_series.reindex(matrix.index).fillna("unknown")
    )

    n_samples, n_features = matrix.shape

    # figure dimensions scale with data size, capped to avoid unreadable plots
    fig_h = max(5, n_samples  * 0.45 + 3.5)
    fig_w = max(9, min(n_features * 0.13 + 4, 30))

    show_col_labels = (max_feature_labels > 0) and (n_features <= max_feature_labels)
    col_fontsize    = max(5, min(8, int(200 / n_features))) if show_col_labels else 8

    g = sns.clustermap(
        matrix,
        method           = linkage,
        metric           = metric,
        cmap             = cmap,
        center           = 0,           # centre colourmap at 0 (correct for scaled data)
        row_colors       = row_colours,
        yticklabels      = True,
        xticklabels      = show_col_labels,
        linewidths       = 0,
        figsize          = (fig_w, fig_h),
        cbar_kws         = {"label": "scaled intensity", "shrink": 0.45},
        dendrogram_ratio = (0.12, 0.15),
        colors_ratio     = 0.025,
    )

    # rotate tick labels
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9
    )
    if show_col_labels:
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=col_fontsize
        )

    g.ax_heatmap.set_xlabel("Features", fontsize=10)
    g.ax_heatmap.set_ylabel("Samples",  fontsize=10)

    # group legend (row colour annotation)
    patches = [
        mpatches.Patch(color=col, label=grp)
        for grp, col in colour_map.items()
    ]
    g.ax_heatmap.legend(
        handles=patches, title="Group",
        bbox_to_anchor=(1.18, 1.0), loc="upper left",
        borderaxespad=0, framealpha=0.9, fontsize=9,
    )

    g.figure.suptitle(
        f"HCA heatmap  ({linkage} linkage, {metric} distance,  "
        f"{n_samples} samples x {n_features} features)",
        y=1.01, fontsize=11, fontweight="bold",
    )

    g.figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(g.figure)
    return g


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 5: HCA ---------------------------------------------------")

    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed.csv")
    groups_path = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    for p in (matrix_path, groups_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} not found - run the preceding pipeline steps first."
            )

    matrix       = pd.read_csv(matrix_path, index_col="sample")
    group_series = pd.read_csv(groups_path,  index_col="sample")["group"]
    group_series = group_series.reindex(matrix.index).fillna("unknown")

    print(f"  samples x features : {matrix.shape}")
    print(f"  linkage            : {cfg.HCA_LINKAGE}")
    print(f"  metric             : {cfg.HCA_METRIC}")
    print(f"  colormap           : {cfg.HCA_CMAP}")

    # build display matrix (compound names when FEATURE_LABEL = "name")
    label_map = _load_feature_labels(cfg)
    if label_map:
        display_matrix = matrix.copy()
        display_matrix.columns = [label_map.get(fid, fid)
                                   for fid in matrix.columns]
        print(f"  feature labels     : compound names  (FEATURE_LABEL='name')")
    else:
        display_matrix = matrix

    # --- clustered heatmap ----------------------------------------------------
    out_plot = os.path.join(plots_dir, "hca_heatmap.png")
    g = plot_hca(
        display_matrix,
        group_series,
        linkage            = cfg.HCA_LINKAGE,
        metric             = cfg.HCA_METRIC,
        cmap               = cfg.HCA_CMAP,
        max_feature_labels = cfg.HCA_MAX_FEATURE_LABELS,
        output_path        = out_plot,
    )
    print(f"  -> {out_plot}")

    # --- save dendrogram orders -----------------------------------------------
    row_idx   = g.dendrogram_row.reordered_ind
    col_idx   = g.dendrogram_col.reordered_ind

    row_order = pd.DataFrame({
        "sample":    matrix.index[row_idx],
        "hca_order": range(len(row_idx)),
    }).set_index("sample")
    out_row = os.path.join(cfg.OUTPUT_DIR, "hca_row_order.csv")
    row_order.to_csv(out_row)
    print(f"  -> {out_row}  (sample dendrogram order)")

    col_order = pd.DataFrame({
        "feature_id": matrix.columns[col_idx],
        "hca_order":  range(len(col_idx)),
    }).set_index("feature_id")
    out_col = os.path.join(cfg.OUTPUT_DIR, "hca_col_order.csv")
    col_order.to_csv(out_col)
    print(f"  -> {out_col}  (feature dendrogram order)")

    return g, row_order, col_order


if __name__ == "__main__":
    run()
