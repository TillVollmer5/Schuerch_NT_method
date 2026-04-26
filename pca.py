"""
pca.py - Step 4 of the GCMS processing pipeline.

Performs Principal Component Analysis on the processed peak matrix and writes:

  output/pca_scores.csv       - sample scores for all N_COMPONENTS
  output/pca_loadings.csv     - feature loadings for all N_COMPONENTS
  output/pca_variance.csv     - explained variance per component
  output/plots/pca_scores.png - scores scatter plot (PC_x vs PC_y)
  output/plots/pca_loadings.png - loadings scatter plot (PC_x vs PC_y)

Axis labels carry the percentage of explained variance, e.g. "PC1 (42.3 %)".
Confidence ellipses (95 %) are drawn per group when PCA_ELLIPSE = True.
The top PCA_TOP_LOADINGS features (by distance from origin) are labelled
in the loadings plot.

The number of components, which PCs to plot, and annotation options are all
controlled through config.py.

Input  : output/peak_matrix_processed.csv
         output/sample_groups.csv
Output : see above

Usage:
    python pca.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive backend - no display needed
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sklearn.decomposition import PCA

import config


# --- Class annotation helpers -------------------------------------------------

def _load_class_annotation(cfg, feature_ids):
    """
    Build two dicts from feature_metadata_enriched.csv based on config:

    highlight_map : feature_id -> hex colour string
        Populated from CLASS_HIGHLIGHT entries (list of dicts with keys
        'column', 'value', 'color').  Features matching multiple entries use
        the last matching color.  Empty dict when CLASS_HIGHLIGHT is empty or
        the enriched file is absent.

    class_label_map : feature_id -> class string
        Populated from CLASS_LABEL_COLUMN.  Only entries with a non-NaN value
        are included.  Empty dict when CLASS_LABEL_COLUMN is "" or file absent.

    Parameters
    ----------
    cfg         : config module
    feature_ids : iterable  all feature_ids present in the loadings DataFrame
    """
    highlight_rules = getattr(cfg, "CLASS_HIGHLIGHT", [])
    label_col       = getattr(cfg, "CLASS_LABEL_COLUMN", "")
    enriched_path   = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")

    highlight_map   = {}
    class_label_map = {}

    if not highlight_rules and not label_col:
        return highlight_map, class_label_map

    if not os.path.exists(enriched_path):
        print("  [info] feature_metadata_enriched.csv not found; "
              "CLASS_HIGHLIGHT and CLASS_LABEL_COLUMN have no effect.")
        return highlight_map, class_label_map

    enriched = pd.read_csv(enriched_path, index_col="feature_id")

    # build highlight_map
    for rule in highlight_rules:
        col   = rule.get("column", "")
        val   = rule.get("value",  "")
        color = rule.get("color",  "#e74c3c")
        if col not in enriched.columns:
            print(f"  [info] CLASS_HIGHLIGHT column '{col}' not in enriched metadata; skipped.")
            continue
        for fid in feature_ids:
            if fid in enriched.index and enriched.loc[fid, col] == val:
                highlight_map[fid] = color

    # build class_label_map
    if label_col:
        if label_col not in enriched.columns:
            print(f"  [info] CLASS_LABEL_COLUMN '{label_col}' not in enriched metadata; skipped.")
        else:
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


# --- Colour / marker palette (extendable for more groups) --------------------

_PALETTE = [
    ("#2166ac", "o"),   # blue   circle
    ("#d6604d", "s"),   # red    square
    ("#4dac26", "^"),   # green  triangle up
    ("#8073ac", "D"),   # purple diamond
    ("#f4a582", "v"),   # orange triangle down
    ("#1a1a1a", "P"),   # black  plus-filled
]

def _group_styles(groups):
    """Return {group: (colour, marker)} for each unique group."""
    unique = list(dict.fromkeys(groups))   # preserve encounter order
    return {g: _PALETTE[i % len(_PALETTE)] for i, g in enumerate(unique)}


# --- PCA computation ---------------------------------------------------------

def fit_pca(matrix, n_components):
    """
    Fit PCA on *matrix* (samples x features).

    Returns
    -------
    scores_df   : DataFrame  samples x PCs
    loadings_df : DataFrame  features x PCs
    variance_df : DataFrame  PC -> explained_variance, ratio, cumulative_ratio
    """
    n_components = min(n_components, *matrix.shape)
    # random_state=0 ensures reproducible results across devices/platforms
    # (guards against the randomized SVD solver being selected for larger matrices)
    pca          = PCA(n_components=n_components, random_state=0)
    scores_arr   = pca.fit_transform(matrix.values)

    pc_names    = [f"PC{i + 1}" for i in range(n_components)]
    scores_df   = pd.DataFrame(scores_arr,   index=matrix.index,   columns=pc_names)
    loadings_df = pd.DataFrame(pca.components_.T,
                               index=matrix.columns, columns=pc_names)
    loadings_df.index.name = "feature_id"

    variance_df = pd.DataFrame({
        "explained_variance":            pca.explained_variance_,
        "explained_variance_ratio":      pca.explained_variance_ratio_,
        "cumulative_explained_variance": np.cumsum(pca.explained_variance_ratio_),
    }, index=pd.Index(pc_names, name="PC"))

    return scores_df, loadings_df, variance_df


# --- Group separation score ---------------------------------------------------

def _group_separation_scores(loadings_df, scores_df, group_series, pc_cols):
    """
    Score each feature by how well its loading direction separates sample groups
    in the score plot.

    For each feature j, computes the quadratic form  l_j^T * B * l_j  where B
    is the between-group scatter matrix of the PCA scores on *pc_cols* and l_j
    is the feature's loading vector on those axes.  Features whose loading
    direction aligns with the axis of maximum group separation score highest.

    Works for any number of groups and any number of PC axes (2D or 3D).

    Parameters
    ----------
    loadings_df  : DataFrame  features x PCs (must contain pc_cols)
    scores_df    : DataFrame  samples  x PCs (must contain pc_cols)
    group_series : Series     sample -> group label (index aligned with scores_df)
    pc_cols      : list[str]  e.g. ["PC1", "PC2"] or ["PC1", "PC2", "PC3"]

    Returns
    -------
    ndarray  shape (n_features,)  — one separation score per feature row
    """
    n_pcs         = len(pc_cols)
    groups        = group_series.reindex(scores_df.index).fillna("unknown").values
    unique_groups = np.unique(groups)
    n_total       = len(groups)

    grand_means = np.array([scores_df[c].mean() for c in pc_cols])

    # Between-group scatter matrix  B  (n_pcs × n_pcs)
    B = np.zeros((n_pcs, n_pcs))
    for g in unique_groups:
        mask = groups == g
        n_g  = int(mask.sum())
        c    = np.array([scores_df.loc[mask, col].mean() for col in pc_cols]) - grand_means
        B   += n_g * np.outer(c, c)
    B /= max(n_total, 1)

    # Loading matrix  L  (n_features × n_pcs)
    L = loadings_df[pc_cols].values

    # Quadratic form for every feature simultaneously:
    #   score_j = l_j^T B l_j  =  diag(L B L^T)
    LB         = L @ B                    # (n_features × n_pcs)
    sep_scores = (LB * L).sum(axis=1)    # elementwise product then row-sum

    return sep_scores


# --- Confidence ellipse helper ------------------------------------------------

def _confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """
    Draw a covariance confidence ellipse for points (x, y) on *ax*.
    *n_std* controls the radius in standard deviations (~95 % for n_std=2).
    Requires at least 2 points; silently skipped otherwise.
    """
    if len(x) < 2:
        return
    from matplotlib.patches import Ellipse

    cov  = np.cov(x, y)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    angle  = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    width  = 2 * n_std * np.sqrt(eigvals[0])
    height = 2 * n_std * np.sqrt(eigvals[1])

    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width, height=height, angle=angle,
        **kwargs,
    )
    ax.add_patch(ellipse)


# --- Plotting -----------------------------------------------------------------

def plot_scores(scores_df, group_series, variance_df, pc_x, pc_y,
                output_path, draw_ellipse=True):
    """
    Scatter plot of PCA scores coloured and marked by group.

    Parameters
    ----------
    scores_df    : DataFrame  samples x PCs
    group_series : Series     sample -> group label
    variance_df  : DataFrame  from fit_pca
    pc_x, pc_y   : int        1-indexed PC numbers for the axes
    output_path  : str
    draw_ellipse : bool        draw 95 % confidence ellipses per group
    """
    col_x  = f"PC{pc_x}"
    col_y  = f"PC{pc_y}"
    pct_x  = variance_df.loc[col_x, "explained_variance_ratio"] * 100
    pct_y  = variance_df.loc[col_y, "explained_variance_ratio"] * 100
    styles = _group_styles(group_series.values)

    fig, ax = plt.subplots(figsize=(6, 5))

    for group, (colour, marker) in styles.items():
        mask = group_series == group
        xs   = scores_df.loc[mask, col_x].values
        ys   = scores_df.loc[mask, col_y].values

        ax.scatter(xs, ys, color=colour, marker=marker, s=70,
                   label=group, zorder=3, edgecolors="white", linewidths=0.5)

        if draw_ellipse and len(xs) >= 2:
            _confidence_ellipse(xs, ys, ax, n_std=2.0,
                                facecolor=colour, alpha=0.12,
                                edgecolor=colour, linewidth=1.2, linestyle="--")

    # sample name labels
    for sample in scores_df.index:
        ax.annotate(sample,
                    (scores_df.loc[sample, col_x], scores_df.loc[sample, col_y]),
                    textcoords="offset points", xytext=(5, 4),
                    fontsize=6.5, color="#444444")

    ax.axhline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.set_xlabel(f"{col_x}  ({pct_x:.1f} %)", fontsize=11)
    ax.set_ylabel(f"{col_y}  ({pct_y:.1f} %)", fontsize=11)
    ax.set_title("PCA scores", fontsize=12, fontweight="bold")
    ax.legend(title="Group", framealpha=0.9, fontsize=9)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_loadings(loadings_df, variance_df, pc_x, pc_y,
                  output_path, top_n=10,
                  highlight_map=None, class_label_map=None,
                  sep_scores=None):
    """
    Scatter plot of feature loadings.
    The *top_n* features are labelled.  When *sep_scores* is provided the
    selection is by group-separation score (between-group scatter quadratic
    form); otherwise Euclidean distance from the origin is used as fallback.

    Parameters
    ----------
    loadings_df     : DataFrame  features x PCs
    variance_df     : DataFrame  from fit_pca
    pc_x, pc_y      : int        1-indexed PC numbers for the axes
    output_path     : str
    top_n           : int        number of features to label (0 = none)
    highlight_map   : dict or None  feature_id -> hex color; highlighted dots
                                    are drawn in their assigned color instead
                                    of the default grey
    class_label_map : dict or None  feature_id -> class string; when present,
                                    "[class]" is appended to the feature label
    sep_scores      : ndarray or None  per-feature group-separation scores
                                    (same row order as loadings_df)
    """
    col_x  = f"PC{pc_x}"
    col_y  = f"PC{pc_y}"
    pct_x  = variance_df.loc[col_x, "explained_variance_ratio"] * 100
    pct_y  = variance_df.loc[col_y, "explained_variance_ratio"] * 100

    lx  = loadings_df[col_x].values
    ly  = loadings_df[col_y].values
    ids = loadings_df.index.tolist()

    highlight_map   = highlight_map   or {}
    class_label_map = class_label_map or {}

    fig, ax = plt.subplots(figsize=(10, 5))

    # draw background (non-highlighted) dots
    # split into highlighted and non-highlighted for correct z-ordering
    hi_idx = [i for i, fid in enumerate(ids) if fid in highlight_map]
    bg_idx = [i for i, fid in enumerate(ids) if fid not in highlight_map]

    if bg_idx:
        ax.scatter(lx[bg_idx], ly[bg_idx], color="#555555", s=18,
                   alpha=0.55, zorder=3, edgecolors="none")
    for i in hi_idx:
        ax.scatter(lx[i], ly[i], color=highlight_map[ids[i]], s=28,
                   alpha=0.85, zorder=4, edgecolors="none")

   ## build legend entries for highlighted classes
   #if highlight_map:
   #    from matplotlib.lines import Line2D
   #    seen = {}
   #    for fid, col in highlight_map.items():
   #        # use the class value as legend label if available
   #        lbl = class_label_map.get(fid, fid)
   #        # group by color
   #        seen.setdefault(col, lbl)
   #    handles = [Line2D([], [], marker="o", color="w",
   #                      markerfacecolor=col, markersize=7, label=lbl)
   #               for col, lbl in seen.items()]
   #    ax.legend(handles=handles, fontsize=7.5, framealpha=0.9,
   #              title="Highlighted", title_fontsize=8,
   #              loc="upper left", bbox_to_anchor=(1.02, 1),
   #              borderaxespad=0)
        
    # label top features by group-separation score (or Euclidean fallback)
    if top_n > 0:
        dist    = sep_scores if sep_scores is not None else np.sqrt(lx ** 2 + ly ** 2)
        top_idx = np.argsort(dist)[-top_n:]
        top_colors = [highlight_map.get(ids[i], "#555555") for i in top_idx]
        ax.scatter(lx[top_idx], ly[top_idx], color=top_colors,
                   s=38, zorder=5, edgecolors="#111111", linewidths=0.9)

        # sort top features bottom-to-top by y so label order matches point order
        top_idx = top_idx[np.argsort(ly[top_idx])]

        xlim  = ax.get_xlim()
        x_off = (xlim[1] - xlim[0]) * 0.25   # small offset in data coords

        # place labels right of point when x >= 0, left of point when x < 0
        texts = []
        for i in top_idx:
            fid   = loadings_df.index[i]
            label = str(fid)
            if fid in class_label_map:
                label = f"{label} [{class_label_map[fid]}]"
            color = highlight_map.get(fid, "#555555")
            if lx[i] >= 0:
                t = ax.text(lx[i] + x_off, ly[i], label,
                            fontsize=6.5, color="#111111", zorder=5,
                            va="center", ha="left")
            else:
                t = ax.text(lx[i] - x_off, ly[i], label,
                            fontsize=6.5, color="#111111", zorder=5,
                            va="center", ha="right")
            texts.append((t, lx[i], ly[i], color))

        # render once to measure real text heights in data coordinates
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        inv      = ax.transData.inverted()

        items = []
        for t, px, py, color in texts:
            bb         = t.get_window_extent(renderer=renderer)
            corners    = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
            h          = abs(corners[1, 1] - corners[0, 1])
            right_side = (px >= 0)
            items.append({"t": t, "px": px, "py": py, "ty": py, "h": h,
                          "right_side": right_side, "color": color})

        # sort bottom-to-top, then push overlapping labels upward
        items.sort(key=lambda d: d["ty"])
        for k in range(1, len(items)):
            min_y = items[k - 1]["ty"] + items[k - 1]["h"] * 1.2
            if items[k]["ty"] < min_y:
                items[k]["ty"] = min_y

        # apply nudged y positions
        for d in items:
            x_cur, _ = d["t"].get_position()
            d["t"].set_position((x_cur, d["ty"]))

        # expand plot limits so all labels fit without flipping
        fig.canvas.draw()
        new_xmin, new_xmax = ax.get_xlim()
        new_ymin, new_ymax = ax.get_ylim()
        for d in items:
            bb = d["t"].get_window_extent(renderer=renderer)
            corners = inv.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
            new_xmin = min(new_xmin, corners[0, 0])
            new_xmax = max(new_xmax, corners[1, 0])
            new_ymin = min(new_ymin, corners[0, 1])
            new_ymax = max(new_ymax, corners[1, 1])
        margin_x = (new_xmax - new_xmin) * 0.03
        margin_y = (new_ymax - new_ymin) * 0.03
        ax.set_xlim(new_xmin - margin_x, new_xmax + margin_x)
        ax.set_ylim(new_ymin - margin_y, new_ymax + margin_y)

        # draw connector line from point to nearest edge of label bbox
        fig.canvas.draw()
        inv = ax.transData.inverted()   # refresh after axis limits changed
        for d in items:
            bb = d["t"].get_window_extent(renderer=renderer)
            if d["right_side"]:
                lx_l, ly_l = inv.transform((bb.x0, (bb.y0 + bb.y1) / 2))
            else:
                lx_l, ly_l = inv.transform((bb.x1, (bb.y0 + bb.y1) / 2))
            ax.plot([d["px"], lx_l], [d["py"], ly_l],
                    color=d["color"], lw=0.5, zorder=3)

    ax.axhline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.set_xlabel(f"{col_x}  ({pct_x:.1f} %)", fontsize=11)
    ax.set_ylabel(f"{col_y}  ({pct_y:.1f} %)", fontsize=11)
    ax.set_title("PCA loadings", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_aspect("equal", adjustable="datalim")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# --- Loading bar chart -------------------------------------------------------

def plot_loadings_bar(loadings_df, variance_df, pc_x, pc_y,
                      output_path, top_n=10,
                      highlight_map=None, class_label_map=None,
                      sep_scores=None):
    """
    Grouped horizontal bar chart showing the top *top_n* features.
    When *sep_scores* is provided, features are selected by group-separation
    score (between-group scatter quadratic form); otherwise Euclidean distance
    from the origin is used as fallback.

    Each feature gets one bar per PC so that positive and negative impacts
    on both components can be read simultaneously.  Features are sorted by
    their PC_x loading so the direction of effect is immediately visible.

    Parameters
    ----------
    loadings_df     : DataFrame  features x PCs
    variance_df     : DataFrame  from fit_pca
    pc_x, pc_y      : int        1-indexed PC numbers (matching PCA_PLOT_X/Y)
    output_path     : str
    top_n           : int        number of features to show (config: PCA_BAR_TOP)
    highlight_map   : dict or None  feature_id -> hex color; y-tick labels for
                                    highlighted features are drawn in that color
    class_label_map : dict or None  feature_id -> class string; appended as
                                    "[class]" to y-tick labels
    sep_scores      : ndarray or None  per-feature group-separation scores
                                    (same row order as loadings_df)
    """
    col_x = f"PC{pc_x}"
    col_y = f"PC{pc_y}"
    pct_x = variance_df.loc[col_x, "explained_variance_ratio"] * 100
    pct_y = variance_df.loc[col_y, "explained_variance_ratio"] * 100

    lx   = loadings_df[col_x].values
    ly   = loadings_df[col_y].values
    dist = sep_scores if sep_scores is not None else np.sqrt(lx ** 2 + ly ** 2)

    highlight_map   = highlight_map   or {}
    class_label_map = class_label_map or {}

    # select top_n by separation score (or Euclidean fallback), then sort by PC_x loading for display
    top_n    = min(top_n, len(loadings_df))
    top_idx  = np.argsort(dist)[-top_n:]
    top_idx  = top_idx[np.argsort(lx[top_idx])]   # sort ascending by PC_x

    feat_ids    = [loadings_df.index[i] for i in top_idx]
    feat_labels = []
    for fid in feat_ids:
        lbl = str(fid)
        if fid in class_label_map:
            lbl = f"{lbl} [{class_label_map[fid]}]"
        feat_labels.append(lbl)

    vals_x = lx[top_idx]
    vals_y = ly[top_idx]

    # layout: one bar group per feature, two bars per group
    y_pos   = np.arange(top_n)
    height  = 0.3                         # bar thickness
    col_pcx = config.BAR_TOP_COL_PCX  # green - PC_x (from config)
    col_pcy = config.BAR_TOP_COL_PCY  # orange - PC_y (from config)

    fig_h = max(4.0, top_n * 0.5 + 1.2)  # scale figure height with n features
    fig, ax = plt.subplots(figsize=(8, fig_h))

    bars_x = ax.barh(y_pos + height / 2, vals_x, height,
                     color=col_pcx, alpha=0.82, label=f"{col_x}  ({pct_x:.1f} %)",
                     zorder=3)
    bars_y = ax.barh(y_pos - height / 2, vals_y, height,
                     color=col_pcy, alpha=0.82, label=f"{col_y}  ({pct_y:.1f} %)",
                     zorder=3)

    # value labels at bar ends
    for bar in bars_x:
        w = bar.get_width()
        ax.text(w + (0.002 if w >= 0 else -0.002), bar.get_y() + bar.get_height() / 2,
                f"{w:+.3f}", va="center", ha="left" if w >= 0 else "right",
                fontsize=6.5, color=col_pcx)
    for bar in bars_y:
        w = bar.get_width()
        ax.text(w + (0.002 if w >= 0 else -0.002), bar.get_y() + bar.get_height() / 2,
                f"{w:+.3f}", va="center", ha="left" if w >= 0 else "right",
                fontsize=6.5, color=col_pcy)

    ax.axvline(0, color="#888888", linewidth=0.9, zorder=2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(feat_labels, fontsize=8)

    # color y-tick labels for highlighted features
    for tick, fid in zip(ax.get_yticklabels(), feat_ids):
        if fid in highlight_map:
            tick.set_color(highlight_map[fid])

    ax.set_xlabel("Loading value", fontsize=10)
    criterion = "group separation" if sep_scores is not None else "Euclidean distance"
    ax.set_title(
        f"Top {top_n} features by {criterion}\n"
        f"({col_x}/{col_y} loading plane)",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=9, framealpha=0.9, loc="upper left")
    ax.grid(True, axis="x", linestyle=":", linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    ax.margins(x=0.2)  # small horizontal margin to avoid cutting off bars

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# --- Standalone class highlight legend ----------------------------------------

def plot_class_legend(cfg, output_path):
    """
    Save a standalone PNG legend for all CLASS_HIGHLIGHT color assignments.

    Entries are grouped by annotation column (e.g. "subclass",
    "npclassifier_pathway") with a bold subheading per column, followed by
    one colored swatch + label per defined value.  An "Unclassified" entry
    at the bottom shows the default grey used for non-highlighted features.

    Silently returns without writing a file when CLASS_HIGHLIGHT is empty.

    Visual style matches hca_class_legend.png (_save_class_legend in hca.py).
    """
    from collections import OrderedDict

    rules = getattr(cfg, "CLASS_HIGHLIGHT", [])
    if not rules:
        return

    # group rules by column, preserving first-appearance order
    sections = OrderedDict()
    for rule in rules:
        col   = rule.get("column", "")
        val   = rule.get("value",  "")
        color = rule.get("color",  "#cccccc")
        if not col or not val:
            continue
        sections.setdefault(col, {})[val] = color
    # append "Unclassified" sentinel as its own section
    sections["Unclassified"] = {"non-highlighted features": (0.33, 0.33, 0.33)}

    if not sections:
        return

    # layout constants — identical to hca._save_class_legend
    PATCH_H = 0.30
    TITLE_H = 0.38
    PAD     = 0.20
    LEFT    = 0.15
    FIG_W   = 4.2

    total_h = PAD
    for title, mapping in sections.items():
        total_h += TITLE_H + len(mapping) * PATCH_H + PAD
    total_h = max(total_h, 1.5)

    fig, ax = plt.subplots(figsize=(FIG_W, total_h))
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, total_h)
    ax.axis("off")

    y = total_h - PAD   # start from top, step downward

    for section_title, mapping in sections.items():
        y -= TITLE_H
        ax.text(
            LEFT, y + TITLE_H * 0.35,
            section_title,
            fontsize=10, fontweight="bold", va="center",
        )
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


# --- Interactive 3D plots (only when N_COMPONENTS == 3) ----------------------

def plot_scores_3d(scores_df, group_series, variance_df, output_path):
    """
    Interactive 3D scores plot (PC1 x PC2 x PC3) saved as a self-contained
    HTML file.  Hover tooltips show the sample name and group.
    """
    import plotly.graph_objects as go

    pct = [variance_df.loc[f"PC{i}", "explained_variance_ratio"] * 100
           for i in (1, 2, 3)]

    styles = _group_styles(group_series.values)
    traces = []

    for group, (colour, _) in styles.items():
        mask    = group_series == group
        members = scores_df.loc[mask]
        traces.append(go.Scatter3d(
            x=members["PC1"],
            y=members["PC2"],
            z=members["PC3"],
            mode="markers+text",
            name=group,
            text=members.index.tolist(),
            textposition="top center",
            textfont=dict(size=9),
            marker=dict(size=7, color=colour, opacity=0.85,
                        line=dict(width=0.5, color="white")),
            hovertemplate=(
                "<b>%{text}</b><br>"
                f"Group: {group}<br>"
                "PC1: %{x:.3f}<br>PC2: %{y:.3f}<br>PC3: %{z:.3f}"
                "<extra></extra>"
            ),
        ))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text="PCA scores (3D)", font=dict(size=15)),
        scene=dict(
            xaxis_title=f"PC1  ({pct[0]:.1f} %)",
            yaxis_title=f"PC2  ({pct[1]:.1f} %)",
            zaxis_title=f"PC3  ({pct[2]:.1f} %)",
            xaxis=dict(backgroundcolor="#f7f7f7", gridcolor="white"),
            yaxis=dict(backgroundcolor="#f0f0f0", gridcolor="white"),
            zaxis=dict(backgroundcolor="#ebebeb", gridcolor="white"),
        ),
        legend=dict(title="Group", bgcolor="rgba(255,255,255,0.85)"),
        margin=dict(l=0, r=0, t=50, b=0),
    )
    fig.write_html(output_path, include_plotlyjs="cdn")


def plot_loadings_3d(loadings_df, variance_df, output_path, top_n=10,
                     sep_scores=None):
    """
    Interactive 3D loadings plot (PC1 x PC2 x PC3) saved as a self-contained
    HTML file.  All features are shown as grey dots; the top *top_n* features
    are highlighted and labelled.  When *sep_scores* is provided selection is
    by group-separation score; otherwise 3D Euclidean distance is used.
    Hover tooltips show the feature ID and all three loading values.
    """
    import plotly.graph_objects as go

    pct = [variance_df.loc[f"PC{i}", "explained_variance_ratio"] * 100
           for i in (1, 2, 3)]

    lx = loadings_df["PC1"].values
    ly = loadings_df["PC2"].values
    lz = loadings_df["PC3"].values
    ids = loadings_df.index.tolist()

    dist    = sep_scores if sep_scores is not None else np.sqrt(lx ** 2 + ly ** 2 + lz ** 2)
    top_idx = set(np.argsort(dist)[-top_n:]) if top_n > 0 else set()
    bg_idx  = [i for i in range(len(ids)) if i not in top_idx]
    hi_idx  = list(top_idx)

    traces = []

    # background features
    traces.append(go.Scatter3d(
        x=lx[bg_idx], y=ly[bg_idx], z=lz[bg_idx],
        mode="markers",
        name="features",
        text=[ids[i] for i in bg_idx],
        marker=dict(size=3, color="#888888", opacity=0.4),
        hovertemplate=(
            "<b>%{text}</b><br>"
            "PC1: %{x:.4f}<br>PC2: %{y:.4f}<br>PC3: %{z:.4f}"
            "<extra></extra>"
        ),
    ))

    # top-loading features
    if hi_idx:
        traces.append(go.Scatter3d(
            x=lx[hi_idx], y=ly[hi_idx], z=lz[hi_idx],
            mode="markers+text",
            name=f"top {top_n}",
            text=[ids[i] for i in hi_idx],
            textposition="top center",
            textfont=dict(size=9, color="#111111"),
            marker=dict(size=6, color="#111111", opacity=0.9),
            hovertemplate=(
                "<b>%{text}</b><br>"
                "PC1: %{x:.4f}<br>PC2: %{y:.4f}<br>PC3: %{z:.4f}"
                "<extra></extra>"
            ),
        ))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text="PCA loadings (3D)", font=dict(size=15)),
        scene=dict(
            xaxis_title=f"PC1  ({pct[0]:.1f} %)",
            yaxis_title=f"PC2  ({pct[1]:.1f} %)",
            zaxis_title=f"PC3  ({pct[2]:.1f} %)",
            xaxis=dict(backgroundcolor="#f7f7f7", gridcolor="white"),
            yaxis=dict(backgroundcolor="#f0f0f0", gridcolor="white"),
            zaxis=dict(backgroundcolor="#ebebeb", gridcolor="white"),
        ),
        legend=dict(title="Features", bgcolor="rgba(255,255,255,0.85)"),
        margin=dict(l=0, r=0, t=50, b=0),
    )
    fig.write_html(output_path, include_plotlyjs="cdn")


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 4: PCA ---------------------------------------------------")

    # load inputs
    # PCA uses the exclusion-filtered matrix; HCA and volcano use the full matrix
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_pca.csv")
    groups_path = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"{matrix_path} not found - run normalization.py first.")
    if not os.path.exists(groups_path):
        raise FileNotFoundError(f"{groups_path} not found - run data_import.py first.")

    matrix       = pd.read_csv(matrix_path, index_col="sample")
    group_series = pd.read_csv(groups_path,  index_col="sample")["group"]

    # align group labels to matrix row order (handles any ordering difference)
    group_series = group_series.reindex(matrix.index).fillna("unknown")

    n = cfg.N_COMPONENTS
    print(f"  samples x features : {matrix.shape}")
    print(f"  components         : {n}")

    # fit PCA
    scores_df, loadings_df, variance_df = fit_pca(matrix, n)

    # print variance summary
    for pc, row in variance_df.iterrows():
        print(f"  {pc}: {row['explained_variance_ratio']*100:.1f} %  "
              f"(cumulative: {row['cumulative_explained_variance']*100:.1f} %)")

    # save data files (always use feature_id as index in CSVs)
    out_scores   = os.path.join(cfg.OUTPUT_DIR, "pca_scores.csv")
    out_loadings = os.path.join(cfg.OUTPUT_DIR, "pca_loadings.csv")
    out_variance = os.path.join(cfg.OUTPUT_DIR, "pca_variance.csv")

    scores_df.to_csv(out_scores)
    loadings_df.to_csv(out_loadings)
    variance_df.to_csv(out_variance)

    print(f"  -> {out_scores}")
    print(f"  -> {out_loadings}")
    print(f"  -> {out_variance}")

    # build display loadings (compound names when FEATURE_LABEL = "name")
    label_map = _load_feature_labels(cfg)
    if label_map:
        display_loadings = loadings_df.copy()
        display_loadings.index = [label_map.get(fid, fid)
                                   for fid in loadings_df.index]
        print(f"  feature labels     : compound names  (FEATURE_LABEL='name')")
    else:
        display_loadings = loadings_df

    # load class annotation (highlight colors + class label suffixes)
    # Always uses original feature_ids (before label_map renaming) for lookup
    highlight_map, class_label_map = _load_class_annotation(cfg, loadings_df.index)
    # remap keys to display labels when label_map is active
    if label_map:
        highlight_map   = {label_map.get(k, k): v for k, v in highlight_map.items()}
        class_label_map = {label_map.get(k, k): v for k, v in class_label_map.items()}
    if highlight_map:
        print(f"  class highlight    : {len(highlight_map)} feature(s) highlighted")
    if class_label_map:
        print(f"  class labels       : {len(class_label_map)} feature(s) have class suffix")

    # validate requested PC axes
    pc_x, pc_y = cfg.PCA_PLOT_X, cfg.PCA_PLOT_Y
    if max(pc_x, pc_y) > n:
        raise ValueError(
            f"PCA_PLOT_X={pc_x} / PCA_PLOT_Y={pc_y} exceed N_COMPONENTS={n}. "
            "Increase N_COMPONENTS in config.py."
        )

    # 2D scores plot (always produced)
    out_plot_scores = os.path.join(plots_dir, "pca_scores.png")
    plot_scores(scores_df, group_series, variance_df,
                pc_x, pc_y, out_plot_scores,
                draw_ellipse=cfg.PCA_ELLIPSE)
    print(f"  -> {out_plot_scores}")

    # compute group-separation scores for top-feature selection
    pc_plot_cols  = [f"PC{pc_x}", f"PC{pc_y}"]
    sep_scores_2d = _group_separation_scores(
        loadings_df, scores_df, group_series, pc_plot_cols)
    print(f"  top feature selection : group-separation score "
          f"({f'PC{pc_x}'}/{f'PC{pc_y}'} between-group scatter)")

    # 2D loadings plot (always produced)
    out_plot_loadings = os.path.join(plots_dir, "pca_loadings.png")
    plot_loadings(display_loadings, variance_df,
                  pc_x, pc_y, out_plot_loadings,
                  top_n=cfg.PCA_TOP_LOADINGS,
                  highlight_map=highlight_map,
                  class_label_map=class_label_map,
                  sep_scores=sep_scores_2d)
    print(f"  -> {out_plot_loadings}")

    # loading bar chart (always produced)
    out_bar = os.path.join(plots_dir, "pca_loadings_bar.png")
    plot_loadings_bar(display_loadings, variance_df,
                      pc_x, pc_y, out_bar,
                      top_n=cfg.PCA_BAR_TOP,
                      highlight_map=highlight_map,
                      class_label_map=class_label_map,
                      sep_scores=sep_scores_2d)
    print(f"  -> {out_bar}")

    # standalone class highlight legend (only when CLASS_HIGHLIGHT is non-empty)
    if getattr(cfg, "CLASS_HIGHLIGHT", []):
        out_legend = os.path.join(plots_dir, "class_highlight_legend.png")
        plot_class_legend(cfg, out_legend)
        print(f"  -> {out_legend}")

    # 3D interactive plots (only when exactly 3 components are computed)
    if n == 3:
        out_3d_scores = os.path.join(plots_dir, "pca_scores_3d.html")
        plot_scores_3d(scores_df, group_series, variance_df, out_3d_scores)
        print(f"  -> {out_3d_scores}  (interactive)")

        sep_scores_3d = _group_separation_scores(
            loadings_df, scores_df, group_series, ["PC1", "PC2", "PC3"])
        out_3d_loadings = os.path.join(plots_dir, "pca_loadings_3d.html")
        plot_loadings_3d(display_loadings, variance_df, out_3d_loadings,
                         top_n=cfg.PCA_TOP_LOADINGS,
                         sep_scores=sep_scores_3d)
        print(f"  -> {out_3d_loadings}  (interactive)")

    return scores_df, loadings_df, variance_df


if __name__ == "__main__":
    run()
