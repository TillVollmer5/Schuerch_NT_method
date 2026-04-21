"""
classification.py - Generate classification.csv for manual annotation.

For each feature used in HCA (peak_matrix_processed.csv), outputs one row with:
  feature_id, compound_name, class (empty), Area per sample, note (empty)

class (1-4) and note are left blank for manual filling.

Input:
  output/peak_matrix_processed_hca.csv - defines the HCA feature set (column headers)
  output/peak_matrix_blank_corrected.csv - raw blank-corrected area values
  output/feature_name_map.csv        - feature_id -> compound_name

Output:
  output/classification.csv

Usage:
    python classification.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import numpy as np
import pandas as pd
import config


def run(cfg=config):
    processed_path    = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_hca.csv")
    bc_path           = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    name_map_path     = os.path.join(cfg.OUTPUT_DIR, "feature_name_map.csv")
    out_path          = os.path.join(cfg.OUTPUT_DIR, "classification.csv")

    # Feature IDs used for HCA (columns of the processed matrix)
    processed = pd.read_csv(processed_path, index_col=0)
    feature_ids = list(processed.columns)

    # Compound names: feature_id -> name (missing -> empty string)
    name_map = pd.read_csv(name_map_path, index_col="feature_id")["compound_name"]

    # Blank-corrected areas: index = feature_id, columns = sample names
    bc = pd.read_csv(bc_path, index_col=0)
    sample_names = list(bc.columns)   # e.g. ["S-R1", ..., "S6"]

    # Library match scores (SI, HRF, Delta RI) from best peak per feature
    meta_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")
    meta = pd.read_csv(meta_path, index_col="feature_id")

    # Reference compounds for class=1 assignment
    # feature_id format: "{mean_mz:.4f}_{mean_rt:.4f}"
    ref_path = os.path.join("DATA", "references", "references.csv")
    farn_path = os.path.join("DATA", "references", "farn1-46.csv")
    ref_df = pd.read_csv(ref_path)
    farn_df = pd.read_csv(farn_path)
    ref_rts   = ref_df["Retention Time"].values.astype(float)
    ref_mzs   = ref_df["Reference m/z"].values.astype(float)
    ref_names = ref_df["Component Name"].values
    farn_rts   = farn_df["Retention Time"].values.astype(float)
    farn_mzs   = farn_df["Reference m/z"].values.astype(float)
    farn_names = farn_df["Component Name"].values

    def _ref_info(fid):
        """Return (class, note) for a feature_id based on reference matching."""
        parts = fid.split("_")
        feat_mz = float(parts[0])
        feat_rt = float(parts[1])
        rt_hit = np.abs(ref_rts - feat_rt) <= cfg.RT_MARGIN
        mz_hit = np.abs(ref_mzs - feat_mz) <= cfg.MZ_TOLERANCE
        idx = np.where(rt_hit & mz_hit)[0]
        if idx.size == 0:
            return "", ""
        matched_name = ref_names[idx[0]]
        return 1, f"reference {matched_name}"
    
    def _ref_farn(fid):
        """Return (class, note) for a feature_id based on farnesene matching."""
        parts = fid.split("_")
        feat_mz = float(parts[0])
        feat_rt = float(parts[1])
        rt_hit = np.abs(farn_rts - feat_rt) <= cfg.RT_MARGIN
        mz_hit = np.abs(farn_mzs - feat_mz) <= cfg.MZ_TOLERANCE
        idx = np.where(rt_hit & mz_hit)[0]
        if idx.size == 0:
            return "", ""
        matched_name = farn_names[idx[0]]
        return 2, f"farnesene {matched_name}"

    def _score_class(fid, current_cls, current_note):
        """Return (class, note) for class 2/3/4 based on library match scores."""
        if current_cls == 1:
            return current_cls, current_note
        # Check if compound name contains "peak"
        compound_name = name_map.get(fid, "")
        if "peak" in compound_name.lower():
            return 4, ""
        if fid not in meta.index:
            return 4, ""
        row = meta.loc[fid]
        si       = float(row.get("si_highest_score", 0) or 0)
        rsi      = float(row.get("rsi_highest_score", 0) or 0)
        hrf      = float(row.get("hrf_highest_score", 0) or 0)
        rhrf     = float(row.get("rhrf_highest_score", 0) or 0)
        delta_ri_raw = row.get("delta_ri_highest_score", None)
        try:
            delta_ri = float(delta_ri_raw) if delta_ri_raw not in (None, "", "N/A") else None
        except (ValueError, TypeError):
            delta_ri = None
        #if si >= 800 and hrf >= 90 and delta_ri is not None and abs(delta_ri) <= 50:
        # Level 2: RI matches within 50 units, RHRF > 75, SI > 500, RSI > 600
        if delta_ri is not None and abs(delta_ri) <= 50 and rhrf > 75 and si > 500 and rsi > 600:
            return 2, ""
        #if si >= 700 and hrf >= 80 and delta_ri is not None and abs(delta_ri) <= 100:
        if rhrf > 75 and si > 500 and rsi > 600:
            return 3, ""
        return 4, ""

    rows = []
    for fid in feature_ids:
        compound_name = name_map.get(fid, "")
        if fid in bc.index:
            area_vals = {f"Area_{s}": bc.at[fid, s] for s in sample_names}
        else:
            area_vals = {f"Area_{s}": "" for s in sample_names}
        cls, note = _ref_info(fid)
        if cls == "":
            cls, note = _ref_farn(fid)
        if cls == "":
            cls, note = _score_class(fid, cls, note)
        row = {"feature_id": fid, "compound_name": compound_name, "class": cls, **area_vals, "note": note}
        rows.append(row)

    col_order = ["feature_id", "compound_name", "class"] + [f"Area_{s}" for s in sample_names] + ["note"]
    out = pd.DataFrame(rows, columns=col_order)
    out.to_csv(out_path, index=False)
    print(f"[classification] Wrote {len(rows)} features to {out_path}")

    # --- Add class column to top_features_analysis.csv -------------------------
    top_features_path = os.path.join(cfg.OUTPUT_DIR, "top_features_analysis.csv")
    if os.path.exists(top_features_path):
        top_df = pd.read_csv(top_features_path)
        # Create mapping from feature_id to class
        class_map = {row["feature_id"]: row["class"] for _, row in out.iterrows()}
        # Add class column
        top_df["class"] = top_df["feature_id"].map(class_map)
        # Write back
        top_df.to_csv(top_features_path, index=False)
        print(f"[classification] Added 'class' column to {top_features_path}")


if __name__ == "__main__":
    run(config)