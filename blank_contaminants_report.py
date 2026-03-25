"""
blank_contaminants_report.py - Report features removed by blank correction.

Combines features_removed_blank.csv, feature_name_map.csv, and
feature_peak_log.csv into a single summary table.

Output columns:
  feature_id        - feature identifier
  compound_name     - compound name from feature_name_map.csv (empty if unknown)
  mean_rt           - mean retention time of the feature (3 decimal places)
  mean_mz           - mean m/z of the feature (3 decimal places)
  max_blank_area    - highest area observed across blank files (3 decimal places)
  samples           - comma-separated list of samples where the feature was
                      detected (only peaks whose ref_mz and rt_aligned are
                      within tolerance of the feature's mean values)

Output file: output/blank_contaminants_report.csv

Usage:
    python blank_contaminants_report.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd

import config


def run(cfg=config):
    output_dir = cfg.OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    removed_path  = os.path.join(output_dir, "features_removed_blank.csv")
    name_map_path = os.path.join(output_dir, "feature_name_map.csv")
    peak_log_path = os.path.join(output_dir, "feature_peak_log.csv")
    out_path      = os.path.join(output_dir, "blank_contaminants_report.csv")

    # --- load features removed by blank correction ----------------------------
    removed = pd.read_csv(removed_path, index_col="feature_id")

    if removed.empty:
        print("  [blank_contaminants_report] No features removed by blank correction.")
        removed.reset_index()[["feature_id"]].to_csv(out_path, index=False)
        return

    # round RT, MZ, max_blank_area to 3 decimal places
    removed["mean_rt"]        = removed["mean_rt"].round(3)
    removed["mean_mz"]        = removed["mean_mz"].round(3)
    removed["max_blank_area"] = removed["max_blank_area"].round(3)

    # --- compound names -------------------------------------------------------
    if os.path.exists(name_map_path):
        name_map = pd.read_csv(name_map_path, index_col="feature_id")["compound_name"]
        removed["compound_name"] = removed.index.map(name_map).fillna("")
    else:
        print("  [blank_contaminants_report] feature_name_map.csv not found; names left blank.")
        removed["compound_name"] = ""

    # --- sample list from peak log (RT + MZ verified) -------------------------
    rt_tol = getattr(cfg, "RT_MARGIN",    0.05)
    mz_tol = getattr(cfg, "MZ_TOLERANCE", 0.005)

    if os.path.exists(peak_log_path):
        peak_log = pd.read_csv(peak_log_path)
        peak_log = peak_log[peak_log["area"] > 0]

        # join feature mean_rt / mean_mz onto each peak log row
        peak_log = peak_log.join(
            removed[["mean_rt", "mean_mz"]].rename(
                columns={"mean_rt": "feat_rt", "mean_mz": "feat_mz"}
            ),
            on="feature_id",
        )

        # keep only rows where both RT and MZ are within tolerance
        rt_ok = (peak_log["rt_aligned"] - peak_log["feat_rt"]).abs() <= rt_tol
        mz_ok = (peak_log["ref_mz"]     - peak_log["feat_mz"]).abs() <= mz_tol
        confirmed = peak_log[rt_ok & mz_ok]

        sample_lists = (
            confirmed.groupby("feature_id")["sample"]
            .apply(lambda s: ", ".join(sorted(s.unique())))
        )
        removed["samples"] = removed.index.map(sample_lists).fillna("")
    else:
        print("  [blank_contaminants_report] feature_peak_log.csv not found; samples left blank.")
        removed["samples"] = ""

    # --- assemble output columns in a clean order ----------------------------
    out_df = removed.reset_index()[
        ["feature_id", "compound_name", "mean_rt", "mean_mz", "max_blank_area", "samples"]
    ]

    out_df.to_csv(out_path, index=False)
    print(f"  -> {out_path}  ({len(out_df)} blank-removed features)")


if __name__ == "__main__":
    run()
