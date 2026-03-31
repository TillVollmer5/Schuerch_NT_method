"""
compound_class_plots.py - Step 2d of the GCMS processing pipeline.

Generates pie charts showing the distribution of compound classes across
sample groups, using the enriched metadata produced by Step 2c
(compound_classification.py).

One figure is produced per entry in CLASS_PIE_COLUMNS, per group
(or combined, depending on CLASS_PIE_GROUPS).  Small slices are merged
into "Other" according to CLASS_PIE_MIN_FRACTION.

Input  : output/feature_metadata_enriched.csv
         output/peak_matrix_blank_corrected.csv   (for group-specific detection)
         output/sample_groups.csv
Output : output/plots/class_pie_<column>_<group>.png  (one per column x group)

Usage:
    python compound_class_plots.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import config


# ---------------------------------------------------------------------------
# Colour palette for pie slices (cycles if more slices than colours)
# ---------------------------------------------------------------------------

_PIE_COLORS = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
    "#aecde8", "#ffbe7d", "#d37295", "#a0cbe8", "#fabfd2",
    "#8cd17d", "#b6992d", "#499894", "#e15759", "#86bcb6",
]


def _build_counts(enriched, col, feature_ids, label):
    """
    Count features per class value for a given set of feature_ids.

    Parameters
    ----------
    enriched    : DataFrame  feature_id-indexed enriched metadata
    col         : str        column name to group by
    feature_ids : list       feature IDs to include (group-specific or all)
    label       : str        group label (for printing only)

    Returns
    -------
    Series: class_value -> count, sorted descending, NaN replaced with
    "Unclassified".
    """
    sub = enriched.loc[enriched.index.isin(feature_ids), col].copy()
    sub = sub.fillna("Unclassified")
    counts = sub.value_counts()
    return counts


def _merge_small(counts, min_fraction):
    """
    Merge slices below min_fraction of total into "Other".
    Returns a new Series (may include "Other").
    """
    if min_fraction <= 0.0:
        return counts
    total  = counts.sum()
    keep   = counts[counts / total >= min_fraction]
    other  = counts[counts / total <  min_fraction].sum()
    if other > 0:
        keep = pd.concat([keep, pd.Series({"Other": other})])
    return keep


def _plot_pie(counts, title, output_path):
    """Draw and save a single pie chart from a counts Series."""
    labels  = counts.index.tolist()
    sizes   = counts.values
    colors  = [_PIE_COLORS[i % len(_PIE_COLORS)] for i in range(len(labels))]

    fig, ax = plt.subplots(figsize=(7, 5.5))

    wedges, texts, autotexts = ax.pie(
        sizes,
        labels=None,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p >= 3 else "",
        startangle=140,
        pctdistance=0.78,
        wedgeprops=dict(linewidth=0.6, edgecolor="white"),
    )

    for at in autotexts:
        at.set_fontsize(7.5)

    # legend with counts
    legend_labels = [f"{lbl}  (n={cnt})" for lbl, cnt in zip(labels, sizes)]
    ax.legend(
        wedges, legend_labels,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8,
        framealpha=0.9,
        title="Class",
        title_fontsize=8,
    )

    ax.set_title(title, fontsize=11, fontweight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run(cfg=config):
    """Generate compound class pie charts per column and group."""

    plots_dir   = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 2d: compound class plots -----------------------------------")

    enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
    matrix_path   = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    groups_path   = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    if not os.path.exists(enriched_path):
        print("  [class plots] feature_metadata_enriched.csv not found; "
              "run compound_classification.py first.  Skipping.")
        return

    enriched = pd.read_csv(enriched_path, index_col="feature_id")

    # --- Resolve which columns to plot ---------------------------------------
    pie_columns   = getattr(cfg, "CLASS_PIE_COLUMNS",   ["superclass"])
    min_fraction  = getattr(cfg, "CLASS_PIE_MIN_FRACTION", 0.02)
    detected_only = getattr(cfg, "CLASS_PIE_DETECTED_ONLY", True)
    pie_groups    = getattr(cfg, "CLASS_PIE_GROUPS", "separate")

    valid_columns = [c for c in pie_columns if c in enriched.columns]
    if not valid_columns:
        print(f"  [class plots] None of CLASS_PIE_COLUMNS {pie_columns} found "
              f"in enriched metadata. Available: {list(enriched.columns)}")
        return

    skipped = [c for c in pie_columns if c not in enriched.columns]
    if skipped:
        print(f"  [class plots] Columns not found, skipped: {skipped}")

    # --- Resolve group membership --------------------------------------------
    # group_map: group_name -> list of feature_ids detected in that group
    if detected_only and os.path.exists(matrix_path) and os.path.exists(groups_path):
        matrix       = pd.read_csv(matrix_path, index_col="feature_id")
        group_series = pd.read_csv(groups_path, index_col="sample")["group"]

        available_groups = sorted(group_series.unique().tolist())

        if pie_groups == "separate":
            groups_to_plot = available_groups
        elif pie_groups == "combined":
            groups_to_plot = ["all samples"]
        elif isinstance(pie_groups, list):
            groups_to_plot = [g for g in pie_groups if g in available_groups]
            missing = [g for g in pie_groups if g not in available_groups]
            if missing:
                print(f"  [class plots] Groups not found in data, skipped: {missing}")
        else:
            groups_to_plot = available_groups

        # build feature lists per group
        group_features = {}
        for grp in groups_to_plot:
            if grp == "all samples":
                samples = group_series.index.tolist()
            else:
                samples = group_series[group_series == grp].index.tolist()
            samples = [s for s in samples if s in matrix.columns]
            # features with area > 0 in at least one sample of this group
            detected = (matrix[samples] > 0).any(axis=1)
            group_features[grp] = detected[detected].index.tolist()
    else:
        # no matrix available or detected_only=False: use all features
        all_features    = enriched.index.tolist()
        if pie_groups == "combined" or not os.path.exists(groups_path):
            groups_to_plot  = ["all samples"]
            group_features  = {"all samples": all_features}
        else:
            group_series    = pd.read_csv(groups_path, index_col="sample")["group"]
            available_groups = sorted(group_series.unique().tolist())
            if isinstance(pie_groups, list):
                groups_to_plot = [g for g in pie_groups if g in available_groups]
            elif pie_groups == "separate":
                groups_to_plot = available_groups
            else:
                groups_to_plot = ["all samples"]
            group_features = {g: all_features for g in groups_to_plot}

    print(f"   columns       : {valid_columns}")
    print(f"   groups        : {list(group_features.keys())}")
    print(f"   detected_only : {detected_only}")
    print(f"   min_fraction  : {min_fraction}")

    # --- Generate plots -------------------------------------------------------
    n_plots = 0
    for col in valid_columns:
        for grp, feat_ids in group_features.items():
            counts = _build_counts(enriched, col, feat_ids, grp)
            if counts.empty:
                print(f"   [skip] {col} / {grp}: no data")
                continue

            counts = _merge_small(counts, min_fraction)

            # build safe filename
            safe_col = col.replace(" ", "_").replace("/", "_")
            safe_grp = str(grp).replace(" ", "_").replace("/", "_")
            fname    = f"class_pie_{safe_col}_{safe_grp}.png"
            out_path = os.path.join(plots_dir, fname)

            n_features = sum(len(v) for v in [feat_ids])
            title = (f"Compound class distribution — {col}\n"
                     f"Group: {grp}  (n={len(feat_ids)} features detected)")

            _plot_pie(counts, title, out_path)
            print(f"   -> {out_path}")
            n_plots += 1

    print(f"-- Step 2d complete  ({n_plots} plot(s) written) ------------------")


def main():
    run(config)


if __name__ == "__main__":
    main()
