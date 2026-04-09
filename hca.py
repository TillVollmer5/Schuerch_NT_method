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
from matplotlib.legend import Legend
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


# --- Class annotation helpers -------------------------------------------------

def _load_class_annotation(cfg, feature_ids):
    """
    Build highlight_map and class_label_map from CLASS_HIGHLIGHT / CLASS_LABEL_COLUMN.

    highlight_map    : feature_id -> hex color  (features to color in tick labels)
    class_label_map  : feature_id -> class string  (suffix to append to labels)
    """
    highlight_map   = {}
    class_label_map = {}

    enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
    highlights = getattr(cfg, "CLASS_HIGHLIGHT", [])
    label_col  = getattr(cfg, "CLASS_LABEL_COLUMN", "")

    if not (highlights or label_col) or not os.path.exists(enriched_path):
        return highlight_map, class_label_map

    enriched = pd.read_csv(enriched_path, index_col="feature_id")

    for entry in highlights:
        col   = entry.get("column")
        val   = entry.get("value")
        color = entry.get("color")
        if not (col and val and color):
            continue
        if col not in enriched.columns:
            continue
        for fid in feature_ids:
            if fid in enriched.index and enriched.at[fid, col] == val:
                highlight_map[fid] = color

    if label_col and label_col in enriched.columns:
        for fid in feature_ids:
            if fid in enriched.index:
                v = enriched.at[fid, label_col]
                if pd.notna(v) and str(v).strip() and str(v) != "nan":
                    class_label_map[fid] = str(v).strip()

    return highlight_map, class_label_map


_TAB20_HEX = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
    "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
]

_GRAY_SPECIAL = {"Unknown": "#cccccc", "Unclassified": "#aaaaaa", "Other": "#bbbbbb"}


def _hex_to_rgb01(hex_color):
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _build_col_colors(cfg, feature_ids):
    """
    Build a col_colors DataFrame for seaborn clustermap from enriched metadata.

    Each column in HCA_CLASS_ANNOTATION_COLUMNS becomes one colored annotation
    strip alongside the feature dendrogram.  Colors are resolved from
    CLASS_COLORS first; values not in that dict are auto-assigned from a tab20
    palette in sorted alphabetical order (so the same value always gets the same
    color across plots).

    Returns
    -------
    col_colors_df : pd.DataFrame  (index = feature_ids, cols = annotation cols)
                    or None when HCA_CLASS_ANNOTATION_COLUMNS is empty / file missing
    legend_data   : dict  {col_name: {value: rgb_tuple, ...}}
    """
    ann_cols = getattr(cfg, "HCA_CLASS_ANNOTATION_COLUMNS", [])
    if not ann_cols:
        return None, {}

    enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
    if not os.path.exists(enriched_path):
        print("  [info] feature_metadata_enriched.csv not found; "
              "skipping class annotation strips")
        return None, {}

    enriched     = pd.read_csv(enriched_path, index_col="feature_id")
    class_colors = getattr(cfg, "CLASS_COLORS", {})

    result_cols = {}
    legend_data = {}

    for col in ann_cols:
        if col not in enriched.columns:
            print(f"  [info] HCA annotation column '{col}' not in enriched "
                  "metadata; skipping")
            continue

        vals = (enriched[col]
                .reindex(feature_ids)
                .fillna("Unknown")
                .astype(str)
                .replace("nan", "Unknown"))

        unique_known = sorted(v for v in vals.unique() if v != "Unknown")

        # seed color_map from CLASS_HIGHLIGHT rules that target this column
        # (last rule wins, matching the PCA/volcano behavior)
        highlight_rules = getattr(cfg, "CLASS_HIGHLIGHT", [])
        color_map = {}
        for rule in highlight_rules:
            if rule.get("column") == col and rule.get("value") in unique_known:
                hex_col = rule["color"].lstrip("#")
                rgb = tuple(int(hex_col[i:i+2], 16) / 255.0 for i in (0, 2, 4))
                color_map[rule["value"]] = rgb

        # assign tab20 colors for any values not covered by CLASS_HIGHLIGHT
        unassigned = [v for v in unique_known if v not in color_map]
        palette    = sns.color_palette("tab20", max(len(unassigned), 1))
        for i, v in enumerate(unassigned):
            color_map[v] = tuple(palette[i])

        color_map["Unknown"] = gray

        result_cols[col] = vals.map(color_map)
        legend_data[col] = {v: color_map[v] for v in unique_known if v in color_map}

    if not result_cols:
        return None, {}

    return pd.DataFrame(result_cols, index=feature_ids), legend_data


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


# --- Standalone legend key ----------------------------------------------------

def _save_class_legend(group_colour_map, col_colors_legend, output_path):
    """
    Write a standalone PNG file that maps every color to its class label.

    One section per annotation column (superclass, npclassifier_pathway, …)
    plus a section for sample groups.  This is always readable regardless of
    the heatmap's figure layout.

    Parameters
    ----------
    group_colour_map  : dict  {group_name: color}
    col_colors_legend : dict  {col_name: {value: color}}  (may be empty)
    output_path       : str
    """
    # collect all sections: sample groups + each annotation column
    sections = {}
    if group_colour_map:
        sections["Sample groups"] = group_colour_map
    if col_colors_legend:
        for col_name, val_colors in col_colors_legend.items():
            if val_colors:
                sections[col_name] = val_colors

    if not sections:
        return

    # layout: one column of patches per section, stacked vertically
    PATCH_H   = 0.30   # inches per entry row
    TITLE_H   = 0.38   # inches per section title
    PAD       = 0.20   # inches between sections
    LEFT      = 0.15   # left margin
    FIG_W     = 4.2    # fixed width

    total_h = PAD
    for title, mapping in sections.items():
        total_h += TITLE_H + len(mapping) * PATCH_H + PAD
    total_h = max(total_h, 1.5)

    fig, ax = plt.subplots(figsize=(FIG_W, total_h))
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    y = total_h - PAD  # start from top, step downward

    for section_title, mapping in sections.items():
        y -= TITLE_H
        ax.text(
            LEFT, y + TITLE_H * 0.35,
            section_title,
            fontsize=10, fontweight="bold", va="center",
        )
        # thin rule under the title
        ax.plot([LEFT, FIG_W - LEFT], [y + TITLE_H * 0.05, y + TITLE_H * 0.05],
                color="0.75", lw=0.6)

        for label, color in mapping.items():
            y -= PATCH_H
            rect = mpatches.FancyBboxPatch(
                (LEFT, y + PATCH_H * 0.15),
                PATCH_H * 0.65, PATCH_H * 0.65,
                boxstyle="round,pad=0.02",
                facecolor=color, edgecolor="0.4", linewidth=0.5,
            )
            ax.add_patch(rect)
            ax.text(
                LEFT + PATCH_H * 0.85, y + PATCH_H * 0.5,
                str(label),
                fontsize=8.5, va="center",
            )
        y -= PAD

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


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
             max_feature_labels, output_path,
             col_colors_df=None, col_colors_legend=None,
             feature_display_labels=None, highlight_map=None):
    """
    Draw and save a bidirectional clustered heatmap with group row annotation.

    Parameters
    ----------
    matrix                : DataFrame  samples x features  (processed, scaled)
    group_series          : Series     sample -> group label
    linkage               : str        scipy linkage method
    metric                : str        distance metric
    cmap                  : str        matplotlib/seaborn diverging colormap
    max_feature_labels    : int        label column axis when n_features <= this
    output_path           : str
    col_colors_df         : DataFrame or None
                            feature_id-indexed DataFrame for class annotation strips
    col_colors_legend     : dict or None
                            {col_name: {value: color}} for strip legends
    feature_display_labels: list or None
                            display label per feature (same order as matrix.columns);
                            includes compound names and optional [class] suffix
    highlight_map         : dict or None
                            display_label -> hex color for x-tick label coloring

    Returns
    -------
    g : seaborn.matrix.ClusterGrid
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
        center           = 0,
        row_colors       = row_colours,
        col_colors       = col_colors_df,
        yticklabels      = True,
        xticklabels      = show_col_labels,
        linewidths       = 0,
        figsize          = (fig_w, fig_h),
        cbar_kws         = {"label": "scaled intensity", "shrink": 0.45},
        dendrogram_ratio = (0.12, 0.15),
        colors_ratio     = 0.025,
    )

    # rotate sample (row) tick labels
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9
    )

    # feature (column) tick labels: apply display labels + optional highlighting
    if show_col_labels:
        # seaborn tick labels are the original column names (feature IDs) after
        # clustering reorder; replace with display labels in the same reordered order
        existing_ticks = g.ax_heatmap.get_xticklabels()
        if feature_display_labels is not None:
            # build feature_id -> display_label lookup
            fid_to_display = dict(zip(matrix.columns, feature_display_labels))
            new_labels = [fid_to_display.get(t.get_text(), t.get_text())
                          for t in existing_ticks]
        else:
            new_labels = [t.get_text() for t in existing_ticks]

        g.ax_heatmap.set_xticklabels(
            new_labels, rotation=90, fontsize=col_fontsize
        )

        # color tick labels for highlighted features
        if highlight_map:
            for tick in g.ax_heatmap.get_xticklabels():
                lbl = tick.get_text()
                if lbl in highlight_map:
                    tick.set_color(highlight_map[lbl])

    g.ax_heatmap.set_xlabel("Features", fontsize=10)
    g.ax_heatmap.set_ylabel("Samples",  fontsize=10)

    # --- legends --------------------------------------------------------------
    # We use Legend(...) directly instead of ax.legend() so that each call
    # creates an independent artist — ax.legend() replaces the previous legend
    # and can't stack multiple boxes reliably.
    legend_y_top = 1.0   # axes-coordinate anchor, stepping downward

    def _add_legend(ax, handles, title, y_anchor, fontsize=9):
        leg = Legend(
            ax, handles, [h.get_label() for h in handles],
            title=title,
            bbox_to_anchor=(1.18, y_anchor), loc="upper left",
            borderaxespad=0, framealpha=0.9, fontsize=fontsize,
        )
        ax.add_artist(leg)
        # estimate height consumed: title + entries
        n = len(handles)
        return y_anchor - (0.06 + n * (0.052 if fontsize >= 9 else 0.045))

    # group legend (row colour annotation)
    group_patches = [
        mpatches.Patch(color=col, label=grp)
        for grp, col in colour_map.items()
    ]
    legend_y_top = _add_legend(
        g.ax_heatmap, group_patches, "Group", legend_y_top, fontsize=9
    )

    # class annotation strip legends — stacked below group legend
    if col_colors_legend:
        for col_name, val_colors in col_colors_legend.items():
            if not val_colors:
                continue
            ann_patches = [
                mpatches.Patch(color=col, label=val)
                for val, col in val_colors.items()
            ]
            # truncate to avoid the legend running off the bottom of the figure
            MAX_LEGEND_ENTRIES = 18
            if len(ann_patches) > MAX_LEGEND_ENTRIES:
                ann_patches = ann_patches[:MAX_LEGEND_ENTRIES]
                ann_patches.append(
                    mpatches.Patch(color=(0.9, 0.9, 0.9), label="(+ more)")
                )
            legend_y_top = _add_legend(
                g.ax_heatmap, ann_patches, col_name, legend_y_top, fontsize=8
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

    feature_ids = list(matrix.columns)

    # --- class annotation strips (col_colors) ---------------------------------
    col_colors_df, col_colors_legend = _build_col_colors(cfg, feature_ids)
    ann_cols = getattr(cfg, "HCA_CLASS_ANNOTATION_COLUMNS", [])
    if col_colors_df is not None:
        print(f"  class annotations  : {ann_cols}")

    # --- class highlight + label suffix (CLASS_HIGHLIGHT / CLASS_LABEL_COLUMN) -
    highlight_map, class_label_map = _load_class_annotation(cfg, feature_ids)
    if highlight_map:
        print(f"  class highlights   : {len(highlight_map)} features")

    # --- build display labels (name + optional [class] suffix) ----------------
    label_map = _load_feature_labels(cfg)
    label_col = getattr(cfg, "CLASS_LABEL_COLUMN", "")

    feature_display_labels = []
    display_highlight_map  = {}  # keyed by display label for tick coloring

    for fid in feature_ids:
        display = label_map.get(fid, fid) if label_map else fid
        if class_label_map.get(fid):
            display = f"{display} [{class_label_map[fid]}]"
        feature_display_labels.append(display)
        if fid in highlight_map:
            display_highlight_map[display] = highlight_map[fid]

    if label_map:
        print(f"  feature labels     : compound names  (FEATURE_LABEL='name')")
    if label_col:
        print(f"  class label suffix : {label_col}")

    # --- clustered heatmap ----------------------------------------------------
    out_plot = os.path.join(plots_dir, "hca_heatmap.png")
    g = plot_hca(
        matrix,          # keep feature_id columns so col_colors_df index matches
        group_series,
        linkage               = cfg.HCA_LINKAGE,
        metric                = cfg.HCA_METRIC,
        cmap                  = cfg.HCA_CMAP,
        max_feature_labels    = cfg.HCA_MAX_FEATURE_LABELS,
        output_path           = out_plot,
        col_colors_df         = col_colors_df,
        col_colors_legend     = col_colors_legend,
        feature_display_labels= feature_display_labels,
        highlight_map         = display_highlight_map,
    )
    print(f"  -> {out_plot}")

    # --- standalone color key -------------------------------------------------
    # Build group colour_map here (mirrors _group_colours used inside plot_hca)
    unique_groups = list(dict.fromkeys(
        group_series.reindex(matrix.index).fillna("unknown")
    ))
    group_colour_map = {grp: _PALETTE[i % len(_PALETTE)]
                        for i, grp in enumerate(unique_groups)}
    out_legend = os.path.join(plots_dir, "hca_class_legend.png")
    _save_class_legend(group_colour_map, col_colors_legend or {}, out_legend)
    print(f"  -> {out_legend}  (color key for annotation strips)")

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
