"""
second_targeted_boxplots.py - Step 7c of the GCMS processing pipeline.

Generates one boxplot panel per feature_id found in top_features_analysis.csv,
comparing groups S and S-R.  Layout and visual design match targeted_boxplots.py.

Input
-----
  output/top_features_analysis.csv
    Columns: feature_id  compound_name  RT  area  sample  PC1  PC2

  Groups are derived from the 'sample' column:
    S    — samples whose name starts with "S" (but not "S-R")
    S-R  — samples whose name starts with "S-R"

Panel title uses compound_name when it is short (≤ 20 characters),
otherwise falls back to feature_id.

Statistical test: Mann-Whitney U (non-parametric).

Output
------
  output/plots/second_targeted_boxplots.png

Usage:
    python second_targeted_boxplots.py
    (or called via pipeline.py)
"""

import math
import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as _config_module

GROUPS_ORDERED = ["S", "S-R"]
GROUP_COLORS   = {"S": "#2166ac", "S-R": "#d6604d"}
MAX_NAME_LEN   = 20


# --- Statistical test ---------------------------------------------------------

def _signif(p):
    if p < 0.0001: return "****"
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return "ns"


def _mannwhitney(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return 1.0, "?"
    try:
        _, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(p), _signif(float(p))
    except Exception:
        return 1.0, "?"


# --- Group assignment ---------------------------------------------------------

def _group(sample_name):
    s = str(sample_name)
    if s.startswith("S-R"):
        return "S-R"
    if s.startswith("S"):
        return "S"
    return None


# --- Main run function --------------------------------------------------------

def run(cfg=_config_module):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 7c: second targeted boxplots -----------------------------")

    csv_path = os.path.join(cfg.OUTPUT_DIR, "top_features_analysis.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found — run top_features_analysis.py first.")

    df = pd.read_csv(csv_path)
    df["group"] = df["sample"].apply(_group)
    df = df[df["group"].notna()]

    feature_ids = list(df["feature_id"].unique())
    n_features  = len(feature_ids)

    if n_features == 0:
        print("  [info] no features found — no plot generated.")
        return

    nrows = getattr(cfg, "TARGETED_BOXPLOT_ROWS", 3)
    ncols = math.ceil(n_features / nrows)

    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2.2, nrows * 3.2))
    axs = np.array(axs).flatten()

    for i in range(n_features, len(axs)):
        axs[i].set_visible(False)

    for idx, fid in enumerate(feature_ids):
        sub = df[df["feature_id"] == fid]
        ax  = axs[idx]

        compound_name = str(sub["compound_name"].iloc[0]) if "compound_name" in sub.columns else ""
        title = fid #compound_name if len(compound_name) <= MAX_NAME_LEN else

        group_data = []
        for g in GROUPS_ORDERED:
            vals = pd.to_numeric(
                sub.loc[sub["group"] == g, "area"], errors="coerce"
            ).dropna().values
            group_data.append(vals)

        _, p_str = _mannwhitney(group_data[0], group_data[1])

        bp = ax.boxplot(
            [d if len(d) > 0 else [np.nan] for d in group_data],
            tick_labels=GROUPS_ORDERED,
            patch_artist=True,
            widths=0.5,
        )
        for patch, g in zip(bp["boxes"], GROUPS_ORDERED):
            patch.set_facecolor(GROUP_COLORS[g])
            patch.set_alpha(0.75)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)

        # overlay individual data points
        for i, (vals, g) in enumerate(zip(group_data, GROUPS_ORDERED)):
            if len(vals) > 0:
                x = np.random.normal(i + 1, 0.04, size=len(vals))
                ax.scatter(x, vals, color=GROUP_COLORS[g], edgecolors="black",
                           linewidths=0.5, s=18, zorder=3, alpha=0.85)

        # y-axis padding
        non_empty = [v for v in group_data if len(v) > 0]
        if not non_empty:
            continue
        all_vals = np.concatenate(non_empty)
        if len(all_vals) > 0 and not np.all(np.isnan(all_vals)):
            lo, hi = np.nanmin(all_vals), np.nanmax(all_vals)
            pad = (hi - lo) * 0.15 if hi > lo else abs(hi) * 0.1 + 0.1
            ax.set_ylim(lo - pad, hi + pad)

        ax.set_title(f"{title}\n{p_str}", fontsize=8, pad=3)
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(True, alpha=0.3, axis="y", linewidth=0.5)

    fig.supylabel("peak area (raw)", fontsize=9)
    plt.tight_layout()

    out_png = os.path.join(plots_dir, "second_targeted_boxplots.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}  ({n_features} features)")


if __name__ == "__main__":
    run()
