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

from sklearn.decomposition import PCA

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
    pca          = PCA(n_components=n_components)
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
                  output_path, top_n=10):
    """
    Scatter plot of feature loadings.
    The *top_n* features with the largest Euclidean distance from the origin
    on the plotted plane are labelled.

    Parameters
    ----------
    loadings_df : DataFrame  features x PCs
    variance_df : DataFrame  from fit_pca
    pc_x, pc_y  : int        1-indexed PC numbers for the axes
    output_path : str
    top_n       : int        number of features to label (0 = none)
    """
    col_x  = f"PC{pc_x}"
    col_y  = f"PC{pc_y}"
    pct_x  = variance_df.loc[col_x, "explained_variance_ratio"] * 100
    pct_y  = variance_df.loc[col_y, "explained_variance_ratio"] * 100

    lx = loadings_df[col_x].values
    ly = loadings_df[col_y].values

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(lx, ly, color="#555555", s=18, alpha=0.55, zorder=3,
               edgecolors="none")

    # label top features by distance from origin
    if top_n > 0:
        dist    = np.sqrt(lx ** 2 + ly ** 2)
        top_idx = np.argsort(dist)[-top_n:]
        for i in top_idx:
            ax.annotate(loadings_df.index[i],
                        (lx[i], ly[i]),
                        textcoords="offset points", xytext=(4, 3),
                        fontsize=6.5, color="#c0392b", zorder=4)
        ax.scatter(lx[top_idx], ly[top_idx], color="#c0392b",
                   s=28, zorder=4, edgecolors="none")

    ax.axhline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.set_xlabel(f"{col_x}  ({pct_x:.1f} %)", fontsize=11)
    ax.set_ylabel(f"{col_y}  ({pct_y:.1f} %)", fontsize=11)
    ax.set_title("PCA loadings", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# --- Loading bar chart -------------------------------------------------------

def plot_loadings_bar(loadings_df, variance_df, pc_x, pc_y,
                      output_path, top_n=10):
    """
    Grouped horizontal bar chart showing the top *top_n* features by
    Euclidean distance from the origin in the PC_x / PC_y loading plane.

    Each feature gets one bar per PC so that positive and negative impacts
    on both components can be read simultaneously.  Features are sorted by
    their PC_x loading so the direction of effect is immediately visible.

    Selection criterion (Euclidean distance) integrates both PCs: a feature
    that loads strongly on either or both axes will be included.

    Parameters
    ----------
    loadings_df : DataFrame  features x PCs
    variance_df : DataFrame  from fit_pca
    pc_x, pc_y  : int        1-indexed PC numbers (matching PCA_PLOT_X/Y)
    output_path : str
    top_n       : int        number of features to show (config: PCA_BAR_TOP)
    """
    col_x = f"PC{pc_x}"
    col_y = f"PC{pc_y}"
    pct_x = variance_df.loc[col_x, "explained_variance_ratio"] * 100
    pct_y = variance_df.loc[col_y, "explained_variance_ratio"] * 100

    lx   = loadings_df[col_x].values
    ly   = loadings_df[col_y].values
    dist = np.sqrt(lx ** 2 + ly ** 2)

    # select top_n by distance, then sort by PC_x loading for display
    top_n    = min(top_n, len(loadings_df))
    top_idx  = np.argsort(dist)[-top_n:]
    top_idx  = top_idx[np.argsort(lx[top_idx])]   # sort ascending by PC_x

    feat_labels = [loadings_df.index[i] for i in top_idx]
    vals_x      = lx[top_idx]
    vals_y      = ly[top_idx]

    # layout: one bar group per feature, two bars per group
    y_pos   = np.arange(top_n)
    height  = 0.35                         # bar thickness
    col_pcx = "#2166ac"                    # blue  - PC_x
    col_pcy = "#d6604d"                    # red   - PC_y

    fig_h = max(4.0, top_n * 0.55 + 1.2)  # scale figure height with n features
    fig, ax = plt.subplots(figsize=(7, fig_h))

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
    ax.set_xlabel("Loading value", fontsize=10)
    ax.set_title(
        f"Top {top_n} features by loading magnitude\n"
        f"(selected by Euclidean distance in {col_x}/{col_y} space)",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=9, framealpha=0.9, loc="lower right")
    ax.grid(True, axis="x", linestyle=":", linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
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


def plot_loadings_3d(loadings_df, variance_df, output_path, top_n=10):
    """
    Interactive 3D loadings plot (PC1 x PC2 x PC3) saved as a self-contained
    HTML file.  All features are shown as grey dots; the top *top_n* features
    by 3D Euclidean distance from the origin are highlighted and labelled.
    Hover tooltips show the feature ID and all three loading values.
    """
    import plotly.graph_objects as go

    pct = [variance_df.loc[f"PC{i}", "explained_variance_ratio"] * 100
           for i in (1, 2, 3)]

    lx = loadings_df["PC1"].values
    ly = loadings_df["PC2"].values
    lz = loadings_df["PC3"].values
    ids = loadings_df.index.tolist()

    dist    = np.sqrt(lx ** 2 + ly ** 2 + lz ** 2)
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
            textfont=dict(size=9, color="#c0392b"),
            marker=dict(size=6, color="#c0392b", opacity=0.9),
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

    # 2D loadings plot (always produced)
    out_plot_loadings = os.path.join(plots_dir, "pca_loadings.png")
    plot_loadings(display_loadings, variance_df,
                  pc_x, pc_y, out_plot_loadings,
                  top_n=cfg.PCA_TOP_LOADINGS)
    print(f"  -> {out_plot_loadings}")

    # loading bar chart (always produced)
    out_bar = os.path.join(plots_dir, "pca_loadings_bar.png")
    plot_loadings_bar(display_loadings, variance_df,
                      pc_x, pc_y, out_bar,
                      top_n=cfg.PCA_BAR_TOP)
    print(f"  -> {out_bar}")

    # 3D interactive plots (only when exactly 3 components are computed)
    if n == 3:
        out_3d_scores = os.path.join(plots_dir, "pca_scores_3d.html")
        plot_scores_3d(scores_df, group_series, variance_df, out_3d_scores)
        print(f"  -> {out_3d_scores}  (interactive)")

        out_3d_loadings = os.path.join(plots_dir, "pca_loadings_3d.html")
        plot_loadings_3d(display_loadings, variance_df, out_3d_loadings,
                         top_n=cfg.PCA_TOP_LOADINGS)
        print(f"  -> {out_3d_loadings}  (interactive)")

    return scores_df, loadings_df, variance_df


if __name__ == "__main__":
    run()
