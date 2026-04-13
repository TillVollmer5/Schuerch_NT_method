"""
compound_correlation.py - Pairwise Pearson correlation scatter grid.

Reads the targeted compound areas from output/targeted_boxplot.csv and
produces a grid of scatter plots:

  Columns (x-axis panels) : Indole  +  all compounds whose name contains "Hex"
  Rows    (y-axis panels)  : every compound in targeted_boxplot.csv

Each panel shows the 12 sample points (S coloured coral, S-R coloured steelblue),
the Pearson r and p-value, and a linear regression line.

Output : output/plots/compound_correlation.png

Usage:
    python compound_correlation.py
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
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

import config


def _pvalue_label(p):
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def run(cfg=config):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv")
    if not os.path.exists(csv_path):
        print(f"  [compound_correlation] {csv_path} not found; "
              "run targeted_boxplots.py first.")
        return

    df = pd.read_csv(csv_path)

    # --- Identify sample columns and group membership -----------------------
    meta_cols   = {"rt", "mz", "name"}
    sample_cols = [c for c in df.columns if c not in meta_cols]
    sr_cols     = [c for c in sample_cols if c.startswith("S-R")]
    s_cols      = [c for c in sample_cols if c.startswith("S") and not c.startswith("S-R")]

    # set compound name as index; data is compounds x samples
    data = df.set_index("name")[sample_cols]

    all_compounds = data.index.tolist()

    # --- Select x-axis compounds (Indole + Hex-containing names) ------------
    x_compounds = [n for n in all_compounds
                   if "indole" in n.lower() or "hex" in n.lower()]

    if not x_compounds:
        print("  [compound_correlation] No Indole or Hex compounds found in targeted_boxplot.csv.")
        return

    y_compounds = all_compounds   # all compounds on y-axis rows

    n_rows = len(y_compounds)
    n_cols = len(x_compounds)

    print(f"  x-axis compounds ({n_cols}): {x_compounds}")
    print(f"  y-axis compounds ({n_rows}): {y_compounds}")

    # --- Build figure --------------------------------------------------------
    fig, axs = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 2.8, n_rows * 2.5),
        squeeze=False,
    )

    # Group colours
    group_colors = {c: "coral"       for c in s_cols}
    group_colors.update({c: "steelblue" for c in sr_cols})

    for row_idx, y_name in enumerate(y_compounds):
        for col_idx, x_name in enumerate(x_compounds):
            ax = axs[row_idx][col_idx]

            # Scatter points coloured by group
            for sc in sample_cols:
                xi = data.loc[x_name, sc]
                yi = data.loc[y_name, sc]
                ax.scatter(xi, yi, color=group_colors[sc], alpha=0.8, s=28,
                           edgecolors="white", linewidths=0.4)

            # Pearson correlation + regression line — per group
            title_parts = []
            for grp_cols, grp_color, grp_label in [
                (s_cols,  "coral",      "S"),
                (sr_cols, "steelblue", "S-R"),
            ]:
                xg = data.loc[x_name, grp_cols].values.astype(float)
                yg = data.loc[y_name, grp_cols].values.astype(float)
                mask = np.isfinite(xg) & np.isfinite(yg)
                xg, yg = xg[mask], yg[mask]
                if len(xg) > 2:
                    r, p = pearsonr(xg, yg)
                    m, b = np.polyfit(xg, yg, 1)
                    x_line = np.linspace(xg.min(), xg.max(), 100)
                    ax.plot(x_line, m * x_line + b,
                            color=grp_color, linewidth=1.2, linestyle="--", alpha=0.9)
                    title_parts.append(f"{grp_label}: r={r:.2f} {_pvalue_label(p)}")
                else:
                    title_parts.append(f"{grp_label}: n/a")
            ax.set_title("\n".join(title_parts), fontsize=6.5, pad=3)

            # Labels only on outer edges
            if row_idx == n_rows - 1:
                ax.set_xlabel(x_name, fontsize=8, labelpad=4)
            else:
                ax.set_xticklabels([])

            if col_idx == 0:
                ax.set_ylabel(y_name, fontsize=8, labelpad=4)
            else:
                ax.set_yticklabels([])

            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.25, linewidth=0.5)

    # --- Legend (S vs S-R) ---------------------------------------------------
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="coral",
               markersize=7, label="S"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="steelblue",
               markersize=7, label="S-R"),
    ]
    fig.legend(handles=legend_elements, loc="upper right",
               fontsize=8, framealpha=0.9, title="Group", title_fontsize=8)

    fig.suptitle("Compound correlations  –  x: Indole & Hex compounds  |  y: all targeted",
                 fontsize=10, fontweight="bold", y=1.002)

    fig.tight_layout()
    out_path = os.path.join(plots_dir, "compound_correlation.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [compound_correlation] saved → {out_path}")


if __name__ == "__main__":
    run()
