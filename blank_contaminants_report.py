"""
blank_contaminants_report.py - Report blank correction decisions with full traceability.

Reads the audit log produced by blank_correction.py and outputs a human-readable
summary of every comparison that resulted in a feature being removed or a sample
cell being zeroed.  Each row corresponds to one (feature, sample_unit, blank)
comparison decision, so the exact blank file and fold change that triggered each
removal can be traced.

Primary input : output/blank_correction_audit.csv  (written by blank_correction.py)
Secondary     : output/feature_name_map.csv         (compound names, if available)

Output : output/blank_contaminants_report.csv

Output columns
--------------
feature_id        - feature identifier
compound_name     - compound name from feature_name_map.csv (empty if unknown)
mean_rt           - mean aligned RT of the feature
mean_mz           - mean reference m/z of the feature
sample_unit       - sample name / group name / "ALL_MEAN" that was compared
unit_type         - "sample" | "group" | "all_mean"
group             - group the sample_unit belongs to
sample_area       - signal area for this comparison unit
blank_name        - blank file (or "BLANK_MAX" / "BLANK_MEAN") used in comparison
blank_area        - blank area used (after m/z gate, if enabled)
blank_rt          - RT of the matched blank peak
blank_mz          - m/z of the matched blank peak
mz_delta          - |feature_mean_mz - blank_mz|
mz_gate_rejected  - True if this blank was rejected by the m/z proximity gate
fold_change       - sample_area / blank_area for this comparison
comparison_failed - True if this specific comparison was below the threshold
fold_threshold    - FOLD_CHANGE_THRESHOLD applied
decision          - "removed" (feature fully dropped) or "zeroed" (cell set to 0)

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


# Output column order
_OUT_COLS = [
    "feature_id", "compound_name",
    "mean_rt", "mean_mz",
    "sample_unit", "unit_type", "group",
    "sample_area",
    "blank_name", "blank_area", "blank_rt", "blank_mz",
    "mz_delta", "mz_gate_rejected",
    "fold_change", "comparison_failed", "fold_threshold",
    "decision",
]


def run(cfg=config):
    output_dir    = cfg.OUTPUT_DIR
    audit_path    = os.path.join(output_dir, "blank_correction_audit.csv")
    name_map_path = os.path.join(output_dir, "feature_name_map.csv")
    out_path      = os.path.join(output_dir, "blank_contaminants_report.csv")

    os.makedirs(output_dir, exist_ok=True)

    # --- prefer new audit log; fall back to legacy if absent ---
    if not os.path.exists(audit_path):
        print("  [blank_contaminants_report] blank_correction_audit.csv not found - "
              "falling back to legacy mode (features_removed_blank.csv).")
        _run_legacy(cfg, out_path, name_map_path)
        return

    audit = pd.read_csv(audit_path)

    # filter to non-kept decisions only
    removed = audit[audit["decision"].isin(["removed", "zeroed"])].copy()

    if removed.empty:
        print("  [blank_contaminants_report] No features removed or zeroed by blank correction.")
        pd.DataFrame(columns=_OUT_COLS).to_csv(out_path, index=False)
        return

    # --- compound names ---
    if os.path.exists(name_map_path):
        name_map = pd.read_csv(name_map_path, index_col="feature_id")["compound_name"]
        removed["compound_name"] = removed["feature_id"].map(name_map).fillna("")
    else:
        removed["compound_name"] = ""

    # --- round numeric columns for readability ---
    for col, dp in [("mean_rt", 4), ("mean_mz", 4),
                    ("blank_rt", 4), ("blank_mz", 4),
                    ("mz_delta", 5), ("sample_area", 1),
                    ("blank_area", 1), ("fold_change", 4)]:
        if col in removed.columns:
            removed[col] = removed[col].round(dp)

    # --- assemble output ---
    out_df = removed[[c for c in _OUT_COLS if c in removed.columns]]
    out_df.to_csv(out_path, index=False)

    n_features = out_df["feature_id"].nunique() if "feature_id" in out_df.columns else 0
    print(f"  -> {out_path}  "
          f"({len(out_df)} comparison row(s) across {n_features} unique feature(s))")


def _run_legacy(cfg, out_path, name_map_path):
    """
    Legacy fallback: read features_removed_blank.csv and feature_peak_log.csv
    (old pipeline output) and produce the old 6-column report format.
    """
    output_dir    = cfg.OUTPUT_DIR
    removed_path  = os.path.join(output_dir, "features_removed_blank.csv")
    peak_log_path = os.path.join(output_dir, "feature_peak_log.csv")

    if not os.path.exists(removed_path):
        print("  [blank_contaminants_report] features_removed_blank.csv not found - nothing to report.")
        pd.DataFrame().to_csv(out_path, index=False)
        return

    removed = pd.read_csv(removed_path, index_col="feature_id")

    if removed.empty:
        print("  [blank_contaminants_report] No features removed by blank correction.")
        removed.reset_index()[["feature_id"]].to_csv(out_path, index=False)
        return

    for col in ["mean_rt", "mean_mz", "max_blank_area"]:
        if col in removed.columns:
            removed[col] = removed[col].round(3)

    if os.path.exists(name_map_path):
        name_map = pd.read_csv(name_map_path, index_col="feature_id")["compound_name"]
        removed["compound_name"] = removed.index.map(name_map).fillna("")
    else:
        removed["compound_name"] = ""

    rt_tol = getattr(cfg, "RT_MARGIN",    0.05)
    mz_tol = getattr(cfg, "MZ_TOLERANCE", 0.005)

    if os.path.exists(peak_log_path):
        peak_log = pd.read_csv(peak_log_path)
        peak_log = peak_log[peak_log["area"] > 0]
        peak_log = peak_log.join(
            removed[["mean_rt", "mean_mz"]].rename(
                columns={"mean_rt": "feat_rt", "mean_mz": "feat_mz"}
            ),
            on="feature_id",
        )
        rt_ok = (peak_log["rt_aligned"] - peak_log["feat_rt"]).abs() <= rt_tol
        mz_ok = (peak_log["ref_mz"]     - peak_log["feat_mz"]).abs() <= mz_tol
        confirmed = peak_log[rt_ok & mz_ok]
        sample_lists = (
            confirmed.groupby("feature_id")["sample"]
            .apply(lambda s: ", ".join(sorted(s.unique())))
        )
        removed["samples"] = removed.index.map(sample_lists).fillna("")
    else:
        removed["samples"] = ""

    out_df = removed.reset_index()[
        ["feature_id", "compound_name", "mean_rt", "mean_mz", "max_blank_area", "samples"]
    ]
    out_df.to_csv(out_path, index=False)
    print(f"  -> {out_path}  ({len(out_df)} blank-removed features, legacy format)")


if __name__ == "__main__":
    run()
