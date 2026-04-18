"""
targeted_boxplots.py - Step 7b of the GCMS processing pipeline.

Generates one boxplot panel per targeted compound, comparing groups defined in
SAMPLE_GROUPS.  Compounds are taken from TARGETED_LIST in config.py (falls back
to EXCLUSION_LIST when TARGETED_LIST is empty).

Feature matching
----------------
Each TARGETED_LIST entry can be:
  [rt, mz, "name"]
      Match by retention time (±EXCLUSION_RT_MARGIN) AND m/z (±EXCLUSION_MZ_TOLERANCE).

  [rt, mz, "name", "feature_id"]
      Match by exact feature_id string (e.g. "41.0384_5.997").
      rt and mz are still displayed; set to None if unknown.

Normalization applied to the blank-corrected matrix before plotting is
controlled by NORMALIZATION_BOXPLOT, LOG_BASE_BOXPLOT, and SCALING_BOXPLOT
in config.py.  The recommended setting is "none" / "none" / "none" so that
the y-axis shows raw chromatographic peak areas — directional biology is
unambiguous and comparable to literature values.

Statistical test: Mann-Whitney U (non-parametric, appropriate for n < 10).

Output
------
  output/targeted_boxplot_data.csv    - wide table of extracted values
  output/plots/targeted_boxplots.png  - boxplot figure

Usage:
    python targeted_boxplots.py
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
from normalization import (
    sum_normalize, median_normalize, pqn_normalize,
    log_transform,
    pareto_scale, auto_scale, vast_scale, range_scale, level_scale,
    _get_settings,
)


# --- Statistical test ---------------------------------------------------------

def _run_stat_test(a, b, test, all_groups=None):
    """
    Run the configured statistical test between group arrays a and b.
    Returns (p_value, label_string).

    test : str  — value of STAT_TEST_BOXPLOT from config
    all_groups : list of arrays — required only for "kruskal" and "anova"
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if test == "mannwhitney":
        _, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    elif test == "ttest":
        _, p = scipy_stats.ttest_ind(a, b, equal_var=False)
    elif test == "ttest_equal":
        _, p = scipy_stats.ttest_ind(a, b, equal_var=True)
    elif test == "wilcoxon":
        if len(a) != len(b):
            return 1.0, "?"   # wilcoxon requires paired samples of equal length
        _, p = scipy_stats.wilcoxon(a, b)
    elif test == "spearman":
        # group index (0 for a, 1 for b) as the x variable
        x = np.array([0] * len(a) + [1] * len(b))
        y = np.concatenate([a, b])
        _, p = scipy_stats.spearmanr(x, y)
    elif test == "pearson":
        x = np.array([0] * len(a) + [1] * len(b))
        y = np.concatenate([a, b])
        _, p = scipy_stats.pearsonr(x, y)
    elif test == "kruskal":
        groups = all_groups if all_groups else [a, b]
        valid  = [g for g in groups if len(g) > 0]
        if len(valid) < 2:
            return 1.0, "?"
        _, p = scipy_stats.kruskal(*valid)
    elif test == "anova":
        groups = all_groups if all_groups else [a, b]
        valid  = [g for g in groups if len(g) > 0]
        if len(valid) < 2:
            return 1.0, "?"
        _, p = scipy_stats.f_oneway(*valid)
    else:
        return 1.0, "?"

    return float(p), _signif(p)


def _signif(p):
    if p < 0.0001: return "****"
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return "ns"


# --- Feature matching ---------------------------------------------------------

def _find_feature(peak_matrix, entry, rt_margin, mz_tol):
    """
    Return a Series of sample values for the compound described by *entry*.

    entry format:
      [rt, mz, name]              - match by RT + mz proximity
      [rt, mz, name, feature_id]  - match by exact feature_id
      [None, None, name, feature_id] - match by exact feature_id only
    """
    # explicit feature_id match
    if len(entry) >= 4 and entry[3]:
        fid = str(entry[3])
        if fid in peak_matrix.index:
            return peak_matrix.loc[fid]
        print(f"    [warn] feature_id '{fid}' not found for '{entry[2]}'")
        return None

    # RT + mz proximity match — return closest RT match within window
    rt_target = entry[0]
    mz_target = entry[1]
    if rt_target is None:
        return None

    best_fid = None
    best_rt_diff = float("inf")
    for fid in peak_matrix.index:
        parts = fid.split("_")
        try:
            mz_f = float(parts[0])
            rt_f = float(parts[1])
        except (ValueError, IndexError):
            continue
        rt_diff = abs(rt_target - rt_f)
        if rt_diff > rt_margin:
            continue
        if mz_target is not None and abs(mz_target - mz_f) > mz_tol:
            continue
        if rt_diff < best_rt_diff:
            best_rt_diff = rt_diff
            best_fid = fid

    return peak_matrix.loc[best_fid] if best_fid is not None else None


# --- Normalization helper -----------------------------------------------------

def _apply_normalization(matrix, norm, log_base, scaling):
    """
    Apply normalization + log + scaling to a features x samples DataFrame.
    Returns a features x samples DataFrame.
    """
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)

    if norm == "sum":
        matrix = sum_normalize(matrix)
    elif norm == "median":
        matrix = median_normalize(matrix)
    elif norm == "pqn":
        matrix = pqn_normalize(matrix)

    matrix = log_transform(matrix, log_base)

    if scaling == "pareto":
        matrix, _ = pareto_scale(matrix)
    elif scaling == "auto":
        matrix, _ = auto_scale(matrix)
    elif scaling == "vast":
        matrix, _ = vast_scale(matrix)
    elif scaling == "range":
        matrix, _ = range_scale(matrix)
    elif scaling == "level":
        matrix, _ = level_scale(matrix)

    return matrix


# --- y-axis label builder -----------------------------------------------------

def _ylabel(norm, log_base, scaling):
    parts = []
    if norm not in ("none", None, ""):
        parts.append(f"{norm}-norm")
    if log_base not in ("none", None, ""):
        lbl = {2: "log2", 10: "log10", math.e: "ln", "e": "ln",
               "sqrt": "sqrt", "cbrt": "cbrt"}.get(log_base, f"log{log_base}")
        parts.append(lbl)
    if scaling not in ("none", None, ""):
        parts.append(f"{scaling}-scaled")
    return "peak area  [" + " + ".join(parts) + "]" if parts else "peak area (raw)"


# --- Main run function --------------------------------------------------------

def run(cfg=_config_module):
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print("-- Step 7b: targeted boxplots ------------------------------------")

    # --- load blank-corrected matrix (features x samples) --------------------
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"{matrix_path} not found — run blank_correction.py first.")
    peak_matrix = pd.read_csv(matrix_path, index_col="feature_id")

    # --- determine targeted compound list ------------------------------------
    targeted_list = getattr(cfg, "TARGETED_LIST", [])
    if not targeted_list:
        targeted_list = getattr(cfg, "EXCLUSION_LIST", [])
    if not targeted_list:
        print("  [info] TARGETED_LIST and EXCLUSION_LIST are both empty — skipping.")
        return

    # --- per-analysis normalization settings ---------------------------------
    norm, log_base, scaling = _get_settings(cfg, "boxplot")
    stat_test = getattr(cfg, "STAT_TEST_BOXPLOT", "mannwhitney")
    print(f"  normalization : {norm}")
    print(f"  log transform : {log_base}")
    print(f"  scaling       : {scaling}")
    print(f"  stat test     : {stat_test}")

    # --- apply normalization to the full matrix first (sample-wise steps need
    #     all samples present for PQN reference spectrum, etc.)
    matrix_norm = _apply_normalization(peak_matrix.copy(), norm, log_base, scaling)

    # --- determine sample columns per group ----------------------------------
    sample_cols = list(matrix_norm.columns)
    group_cols  = {}   # group_name -> [col, ...]
    for group_name, prefix in cfg.SAMPLE_GROUPS:
        group_cols[group_name] = [
            c for c in sample_cols
            if c.startswith(prefix)
            # exclude columns that belong to a more specific (earlier) group
            and not any(
                c.startswith(other_prefix)
                for _, other_prefix in cfg.SAMPLE_GROUPS
                if other_prefix != prefix and len(other_prefix) > len(prefix)
            )
        ]

    groups_ordered = [g for g, _ in cfg.SAMPLE_GROUPS]

    rt_margin = getattr(cfg, "EXCLUSION_RT_MARGIN", 0.05)
    mz_tol    = getattr(cfg, "EXCLUSION_MZ_TOLERANCE", 0.001)

    # --- extract targeted data -----------------------------------------------
    result_rows = []
    for entry in targeted_list:
        name = entry[2] if len(entry) >= 3 else str(entry)
        values = _find_feature(matrix_norm, entry, rt_margin, mz_tol)
        if values is None:
            print(f"  [warn] '{name}' not found in peak matrix — skipping.")
            row = {"rt": entry[0], "mz": entry[1], "name": name}
            row.update({c: np.nan for c in sample_cols})
        else:
            row = {"rt": entry[0], "mz": entry[1], "name": name}
            row.update(values.to_dict())
        result_rows.append(row)

    targeted_df = pd.DataFrame(result_rows)

    # save extracted data
    out_csv = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_data.csv")
    targeted_df.to_csv(out_csv, index=False)
    print(f"  -> {out_csv}")

    # --- build boxplots -------------------------------------------------------
    n_compounds = len(targeted_df)
    if n_compounds == 0:
        print("  [info] no compounds matched — no plot generated.")
        return

    nrows = getattr(cfg, "TARGETED_BOXPLOT_ROWS", 3)
    ncols = math.ceil(n_compounds / nrows)

    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2.2, nrows * 3.2))
    axs = np.array(axs).flatten()

    # hide unused panels
    for i in range(n_compounds, len(axs)):
        axs[i].set_visible(False)

    # color cycle for groups
    default_colors = [    "#2166ac",  "#d6604d"]
    group_colors = {g: default_colors[i % len(default_colors)]
                    for i, g in enumerate(groups_ordered)}

    ylabel = _ylabel(norm, log_base, scaling)

    for idx, row in targeted_df.iterrows():
        name = row["name"]
        ax   = axs[idx]

        group_data  = []
        group_labels = []
        for g in groups_ordered:
            cols = group_cols.get(g, [])
            vals = pd.to_numeric(row[cols], errors="coerce").dropna().values if cols else np.array([])
            group_data.append(vals)
            group_labels.append(g)

        # significance test
        stat_test = getattr(cfg, "STAT_TEST_BOXPLOT", "mannwhitney")
        p_str = ""
        if len(group_data) >= 2 and len(group_data[0]) > 0 and len(group_data[1]) > 0:
            try:
                _, p_str = _run_stat_test(
                    group_data[0], group_data[1],
                    stat_test,
                    all_groups=group_data,
                )
            except Exception:
                p_str = "?"

        bp = ax.boxplot(
            [d if len(d) > 0 else [np.nan] for d in group_data],
            tick_labels=group_labels,
            patch_artist=True,
            widths=0.5,
        )
        for patch, g in zip(bp["boxes"], groups_ordered):
            patch.set_facecolor(group_colors[g])
            patch.set_alpha(0.75)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)

        # overlay individual data points
        #for i, (vals, g) in enumerate(zip(group_data, groups_ordered)):
        #    if len(vals) > 0:
        #        x = np.random.normal(i + 1, 0.04, size=len(vals))
        #        ax.scatter(x, vals, color=group_colors[g], edgecolors="black",
        #                   linewidths=0.5, s=18, zorder=3, alpha=0.85)

        # y-axis padding
        non_empty = [v for v in group_data if len(v) > 0]
        if not non_empty:
            continue
        all_vals = np.concatenate(non_empty)
        if len(all_vals) > 0 and not np.all(np.isnan(all_vals)):
            lo, hi = np.nanmin(all_vals), np.nanmax(all_vals)
            pad = (hi - lo) * 0.15 if hi > lo else abs(hi) * 0.1 + 0.1
            ax.set_ylim(lo - pad, hi + pad)

        ax.set_title(f"{name}\n{p_str}", fontsize=8, pad=3)
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(True, alpha=0.3, axis="y", linewidth=0.5)

    fig.supylabel(ylabel, fontsize=9)
    plt.tight_layout()

    out_png = os.path.join(plots_dir, "targeted_boxplots.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}  ({n_compounds} compounds)")


if __name__ == "__main__":
    run()
