"""
targeted_boxplots_unnorm.py - Unnormalized targeted boxplots (standalone helper).

Reads peak_matrix_blank_corrected.csv from the current run folder, matches the
compounds in EXCLUSION_LIST by RT and m/z, and produces one boxplot grid showing
raw chromatographic areas.

Output
------
  {OUTPUT_DIR}/targeted_boxplot_unnorm.csv
  {OUTPUT_DIR}/plots/targeted_boxplots_unnorm.png

Usage:
    python targeted_boxplots_unnorm.py
    (or call run(cfg) from pipeline.py)
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

import config as _config_module


def run(cfg=_config_module):
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(
            f"{matrix_path} not found — run blank_correction.py first.")

    peak_matrix     = pd.read_csv(matrix_path)
    targeted_list   = cfg.EXCLUSION_LIST
    sample_cols_all = [c for c in peak_matrix.columns if c != "feature_id"]
    column_names    = ["rt", "mz", "name"] + sample_cols_all

    result_rows = []
    for target in targeted_list:
        mz_target = float(target[1])
        rt_target = float(target[0])
        row_data  = list(target[:3])
        for _, feat_row in peak_matrix.iterrows():
            fid        = feat_row["feature_id"].split("_")
            mz_feature = float(fid[0])
            rt_feature = float(fid[1])
            if (abs(rt_target - rt_feature) <= cfg.EXCLUSION_RT_MARGIN and
                    abs(mz_target  - mz_feature) <= cfg.EXCLUSION_MZ_TOLERANCE):
                row_data = list(target[:3]) + [feat_row[c] for c in sample_cols_all]
                break
        result_rows.append(row_data)

    targeted_df = pd.DataFrame(result_rows, columns=column_names)
    out_csv = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_unnorm.csv")
    targeted_df.to_csv(out_csv, index=False)

    sr_cols = [c for c in sample_cols_all if c.startswith("S-R")]
    s_cols  = [c for c in sample_cols_all
                if c.startswith("S") and not c.startswith("S-R")]

    nrows = 3
    ncols = math.ceil(len(targeted_df) / nrows)
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2, nrows * 3))
    axs = np.array(axs).flatten()
    for i in range(len(targeted_df), nrows * ncols):
        axs[i].set_visible(False)

    fig.supylabel("y-axis: unnormalized area")

    for idx, row in targeted_df.iterrows():
        name = row["name"]
        SR   = np.array(row[sr_cols].values, dtype=float)
        S    = np.array(row[s_cols].values,  dtype=float)

        t_stat, p_value = ttest_ind(S, SR)
        if   p_value < 0.0001: p_signif = "****"
        elif p_value < 0.001:  p_signif = "***"
        elif p_value < 0.01:   p_signif = "**"
        elif p_value < 0.05:   p_signif = "*"
        else:                  p_signif = "ns"

        ax = axs[idx]
        bp = ax.boxplot([S, SR], tick_labels=["S", "SR"], patch_artist=True)
        for patch, color in zip(bp["boxes"], ["coral", "steelblue"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for med in bp["medians"]:
            med.set_color("black")

        all_data = np.concatenate([S, SR])
        if len(all_data):
            lo, hi = np.nanmin(all_data), np.nanmax(all_data)
            pad    = (hi - lo) * 0.1 if hi > lo else abs(hi) * 0.1 + 0.1
            ax.set_ylim(lo - pad, hi + pad)

        ax.set_title(f"{name}\n{p_signif}")
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    out_png = os.path.join(plots_dir, "targeted_boxplots_unnorm.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}")


if __name__ == "__main__":
    run()
