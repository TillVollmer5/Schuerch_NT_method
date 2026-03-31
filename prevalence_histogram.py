"""
prevalence_histogram.py - Prevalence distribution of detected features.

Reads the blank-corrected peak matrix and computes, for each feature, the
fraction of samples in which that feature is detected (area > 0).  A histogram
is produced showing how many features share each prevalence level.

Because the finest achievable prevalence step is 1/N_samples, the bin edges are
automatically aligned to those natural discrete steps, keeping the x-axis
meaningful regardless of how many samples are in the experiment.

Vertical reference lines mark the prevalence thresholds configured for PCA and
HCA so you can immediately see how many features each filter retains or removes.

Input  : output/peak_matrix_blank_corrected.csv
Output : output/plots/prevalence_histogram.png
         output/prevalence_summary.csv

Usage:
    python prevalence_histogram.py
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

import config


def run(cfg):
    """Entry point used by pipeline.py (mirrors the run(cfg) convention)."""
    _run(cfg)


def main():
    _run(config)


def _run(cfg):
    # --- Paths ------------------------------------------------------------------
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    plots_dir   = os.path.join(cfg.OUTPUT_DIR, "plots")
    out_plot    = os.path.join(plots_dir, "prevalence_histogram.png")
    out_csv     = os.path.join(cfg.OUTPUT_DIR, "prevalence_summary.csv")

    os.makedirs(plots_dir, exist_ok=True)

    # --- Load matrix ------------------------------------------------------------
    if not os.path.exists(matrix_path):
        sys.exit(
            f"[prevalence_histogram] Input not found: {matrix_path}\n"
            "Run data_import.py and blank_correction.py first."
        )

    # Rows = features (index = feature_id), columns = sample names
    df = pd.read_csv(matrix_path, index_col=0)

    n_features = df.shape[0]
    n_samples  = df.shape[1]

    if n_samples == 0:
        sys.exit("[prevalence_histogram] Matrix has no sample columns.")

    print(f"[prevalence_histogram] {n_features} features x {n_samples} samples")

    # --- Compute prevalence per feature -----------------------------------------
    # A feature is "detected" in a sample when its area is > 0.
    # prevalence = (number of samples with area > 0) / n_samples
    detected   = (df > 0).sum(axis=1)        # Series: count of non-zero cells per row
    prevalence = detected / n_samples         # float in [0.0, 1.0]

    # --- Bin edges aligned to natural 1/N_samples steps -------------------------
    # Possible values: 0/N, 1/N, 2/N, ..., N/N
    # Bin edges are placed at (k - 0.5)/N so each bar is centred on an exact k/N.
    # np.arange(-0.5, n_samples + 1) -> [-0.5, 0.5, 1.5, ..., N+0.5]  (N+2 values)
    step      = 1.0 / n_samples
    bin_edges = np.arange(-0.5, n_samples + 1) * step   # N+1 bins

    counts, _ = np.histogram(prevalence.values, bins=bin_edges)

    # k_values[i] = number of samples in which features counted by counts[i] were detected
    k_values  = np.arange(0, n_samples + 1)
    centres   = k_values / n_samples

    # --- Summary CSV ------------------------------------------------------------
    summary = pd.DataFrame({
        "n_samples_detected": k_values,
        "prevalence":         centres,
        "n_features":         counts,
    })
    summary.to_csv(out_csv, index=False)
    print(f"[prevalence_histogram] Summary written to {out_csv}")

    # --- Plot -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(max(8, n_samples * 0.65 + 3), 5))

    bar_width = step * 0.75

    bars = ax.bar(
        centres,
        counts,
        width=bar_width,
        color="#4C72B0",
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
        label="Features",
    )

    # Count labels above each non-empty bar
    y_max = max(counts) if max(counts) > 0 else 1
    for bar, cnt in zip(bars, counts):
        if cnt > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height() + y_max * 0.012,
                str(cnt),
                ha="center", va="bottom",
                fontsize=7.5, color="#222222",
            )

    # Threshold reference lines from config.py
    thresholds = [
        ("PCA",     cfg.MIN_PREVALENCE_PCA,     "#D93025", "--"),
        ("HCA",     cfg.MIN_PREVALENCE_HCA,     "#E07B00", ":"),
        ("Volcano", cfg.MIN_PREVALENCE_VOLCANO, "#2E8B57", "-."),
    ]
    legend_entries = []
    for name, thresh, color, ls in thresholds:
        if thresh > 0.0:
            line = ax.axvline(
                thresh - step / 2,      # left edge of the first included bin
                color=color, linestyle=ls, linewidth=1.8, zorder=4,
                label=f"{name} min prevalence ({thresh:.0%})",
            )
            legend_entries.append(line)

    # X-axis: one tick per k/N value; label as "k/N (xx%)"
    ax.set_xticks(centres)
    if n_samples <= 20:
        tick_labels = [
            f"{k}/{n_samples}\n({k/n_samples:.0%})" for k in k_values
        ]
        ax.set_xticklabels(tick_labels, fontsize=8)
    else:
        # Thin out labels to avoid overlap when there are many samples
        tick_step  = max(1, round(n_samples / 10))
        visible    = set(range(0, n_samples + 1, tick_step)) | {n_samples}
        tick_labels = [
            f"{k}/{n_samples}\n({k/n_samples:.0%})" if k in visible else ""
            for k in k_values
        ]
        ax.set_xticklabels(tick_labels, fontsize=8)

    ax.set_xlabel("Prevalence  (fraction of samples with area > 0)", fontsize=11)
    ax.set_ylabel("Number of features (peaks)", fontsize=11)
    ax.set_title(
        f"Feature prevalence distribution  "
        f"({n_features} features, {n_samples} samples)",
        fontsize=13, fontweight="bold", pad=10,
    )

    ax.set_xlim(-step * 0.7, 1.0 + step * 0.7)
    ax.set_ylim(0, y_max * 1.14)
    ax.yaxis.grid(True, linestyle="--", alpha=0.45, zorder=0)
    ax.set_axisbelow(True)

    if legend_entries:
        ax.legend(fontsize=9, loc="upper left", framealpha=0.85)

    plt.tight_layout()
    fig.savefig(out_plot, dpi=200)
    plt.close(fig)
    print(f"[prevalence_histogram] Plot saved to {out_plot}")


if __name__ == "__main__":
    main()
