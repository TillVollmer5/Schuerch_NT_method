"""
blank_correction.py - Step 2 of the GCMS processing pipeline.

Removes features (or individual sample cells) whose signal does not exceed the
blank background by at least FOLD_CHANGE_THRESHOLD.  Features absent from all
blanks are always retained.

Sample-side comparison modes (BLANK_SAMPLE_MODE in config.py)
-------------------------------------------------------------
"mean"       - mean area across all samples vs the blank reference (default).
               Features that fail are removed entirely from the matrix.
"per_sample" - each sample's area is compared individually.  Failing cells are
               zeroed in the matrix rather than removing the whole feature row.
               A feature row is only dropped when every sample cell is zeroed.
"per_group"  - mean area per SAMPLE_GROUPS group is compared.  Failing groups
               have all their sample cells zeroed.  Same all-zero -> drop rule.

Blank reference modes (BLANK_REFERENCE_MODE in config.py)
----------------------------------------------------------
"max"  - highest blank area across all blank files within the RT window (default).
"mean" - mean blank area across all blank files (m/z gate applied per blank before
         averaging; a blank that fails the gate contributes 0 to the mean).
"each" - each blank file is compared independently; a sample unit fails if its
         fold change is below FOLD_CHANGE_THRESHOLD for ANY blank file.

Optional m/z gate (BLANK_USE_MZ in config.py)
---------------------------------------------
When enabled, a blank peak is only accepted as a match if its m/z is within
BLANK_MZ_TOLERANCE Da of the feature's mean_mz.  Blanks that fail the m/z gate
have their area set to 0, treating the feature as "not found in that blank".
Gate details are recorded in the audit log (mz_gate_rejected column).

Inputs
------
  output/peak_matrix_raw.csv        - features x samples
  output/blank_per_feature.csv      - per-(feature, blank-file) blank signal
                                      (written by data_import.py)
  output/blank_features.csv         - fallback when blank_per_feature.csv absent
  output/feature_metadata.csv       - mean_rt, mean_mz for each feature
  output/sample_groups.csv          - sample -> group mapping

Outputs
-------
  output/peak_matrix_blank_corrected.csv   - corrected matrix
  output/blank_correction_audit.csv        - full per-(feature, sample_unit, blank)
                                             comparison log with decisions
  output/features_removed_blank.csv        - backward-compat summary of removed features

Usage:
    python blank_correction.py
"""

import math
import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd

import config


# ---------------------------------------------------------------------------
# Helper: m/z gate
# ---------------------------------------------------------------------------

def _apply_mz_gate(blank_df, feat_mz_map, mz_tol):
    """
    Apply the m/z proximity gate to every row of *blank_df*.

    Parameters
    ----------
    blank_df    : DataFrame with columns feature_id, blank_name,
                  blank_area, blank_rt, blank_mz
    feat_mz_map : dict  feature_id -> mean_mz
    mz_tol      : float  maximum accepted |feature_mz - blank_mz|

    Returns a copy with two new columns added and blank_area zeroed where
    the gate rejects the match:
        mz_delta        - |feature_mean_mz - blank_mz| (NaN if blank_mz is NaN)
        mz_gate_rejected - True where the blank peak was rejected
    """
    df = blank_df.copy()
    feat_mz = df["feature_id"].map(feat_mz_map)

    df["mz_delta"] = (feat_mz - df["blank_mz"]).abs()

    has_area    = df["blank_area"] > 0
    mz_rejected = has_area & (
        df["blank_mz"].isna() | (df["mz_delta"] > mz_tol)
    )
    df["mz_gate_rejected"] = mz_rejected
    df.loc[mz_rejected, "blank_area"] = 0.0
    return df


# ---------------------------------------------------------------------------
# Helper: aggregate blank reference
# ---------------------------------------------------------------------------

def _aggregate_blank_reference(blank_df, mode):
    """
    Collapse the per-file blank table into the reference form required by
    the chosen BLANK_REFERENCE_MODE.

    Parameters
    ----------
    blank_df : DataFrame (already m/z-gated) with columns:
               feature_id, blank_name, blank_area, blank_rt, blank_mz,
               mz_delta, mz_gate_rejected
    mode     : "max" | "mean" | "each"

    Returns
    -------
    DataFrame in the same column format.
    - "each"       : blank_df unchanged (one row per feature × blank_name)
    - "max"/"mean" : one row per feature; blank_name set to "BLANK_MAX" or
                     "BLANK_MEAN"; blank_rt/mz from the highest-area row.
    """
    if mode == "each":
        return blank_df.copy()

    rows = []
    for fid, grp in blank_df.groupby("feature_id", sort=False):
        if mode == "max":
            best = grp.loc[grp["blank_area"].idxmax()]
            ref_area         = float(best["blank_area"])
            ref_rt           = best["blank_rt"]
            ref_mz           = best["blank_mz"]
            ref_mz_delta     = best["mz_delta"]
            ref_mz_rejected  = bool(best["mz_gate_rejected"])
            ref_name         = "BLANK_MAX"
        else:  # "mean"
            ref_area = float(grp["blank_area"].mean())
            # metadata from the blank that had the highest original area
            best = grp.loc[grp["blank_area"].idxmax()]
            ref_rt           = best["blank_rt"]
            ref_mz           = best["blank_mz"]
            ref_mz_delta     = best["mz_delta"]
            ref_mz_rejected  = bool(grp["mz_gate_rejected"].any())
            ref_name         = "BLANK_MEAN"

        rows.append({
            "feature_id":       fid,
            "blank_name":       ref_name,
            "blank_area":       ref_area,
            "blank_rt":         ref_rt,
            "blank_mz":         ref_mz,
            "mz_delta":         ref_mz_delta,
            "mz_gate_rejected": ref_mz_rejected,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Helper: sample units
# ---------------------------------------------------------------------------

def _compute_sample_units(matrix, group_map, mode):
    """
    Build the list of comparison units for the given BLANK_SAMPLE_MODE.

    Each unit is a dict with keys:
        sample_unit      : str   label used in the audit log
        unit_type        : str   "all_mean" | "sample" | "group"
        group            : str   group name or "-"
        area_series      : pd.Series indexed by feature_id
        affected_samples : list[str]  sample columns to zero on failure

    Parameters
    ----------
    matrix    : DataFrame  features x samples
    group_map : dict  sample_name -> group_name
    mode      : "mean" | "per_sample" | "per_group"
    """
    if mode == "mean":
        return [{
            "sample_unit":      "ALL_MEAN",
            "unit_type":        "all_mean",
            "group":            "-",
            "area_series":      matrix.mean(axis=1),
            "affected_samples": list(matrix.columns),
        }]

    if mode == "per_sample":
        return [
            {
                "sample_unit":      col,
                "unit_type":        "sample",
                "group":            group_map.get(col, "unknown"),
                "area_series":      matrix[col],
                "affected_samples": [col],
            }
            for col in matrix.columns
        ]

    if mode == "per_group":
        groups = {}
        for col in matrix.columns:
            g = group_map.get(col, "unknown")
            groups.setdefault(g, []).append(col)
        return [
            {
                "sample_unit":      g,
                "unit_type":        "group",
                "group":            g,
                "area_series":      matrix[cols].mean(axis=1),
                "affected_samples": cols,
            }
            for g, cols in groups.items()
        ]

    raise ValueError(f"Unknown BLANK_SAMPLE_MODE: {mode!r}")


# ---------------------------------------------------------------------------
# Core: build audit and apply corrections
# ---------------------------------------------------------------------------

def _build_and_apply(matrix, blank_ref, sample_units, fold_threshold, sample_mode):
    """
    For every (feature, sample_unit, blank_reference) combination, compute the
    fold change, record the decision, and accumulate matrix corrections.

    Decisions
    ---------
    "kept"    - fold_change >= fold_threshold, or no blank signal (always kept)
    "zeroed"  - per_sample or per_group mode: this cell will be zeroed
    "removed" - all_mean mode: entire feature row will be dropped

    Returns
    -------
    corrected_matrix : DataFrame  with corrections applied
    audit_df         : DataFrame  full audit log
    removed_fids     : list[str]  feature_ids fully removed from the matrix
    """
    # index blank_ref by feature_id for O(1) lookup
    blank_by_fid = {}
    for fid, grp in blank_ref.groupby("feature_id", sort=False):
        blank_by_fid[fid] = grp.to_dict("records")

    audit_rows   = []
    cells_to_zero = []          # list of (feature_id, sample_col)
    fids_to_remove = set()      # features removed entirely (all_mean mode)

    for unit in sample_units:
        su_name          = unit["sample_unit"]
        su_type          = unit["unit_type"]
        su_group         = unit["group"]
        area_series      = unit["area_series"]
        affected_samples = unit["affected_samples"]

        for fid in matrix.index:
            sample_area = float(area_series.get(fid, 0.0)
                                if hasattr(area_series, "get")
                                else area_series.loc[fid]
                                if fid in area_series.index else 0.0)
            if sample_area == 0.0:
                # cell already absent — skip blank comparison
                continue

            blank_rows = blank_by_fid.get(fid, [])
            if not blank_rows:
                continue  # no blank data at all — always kept, no audit row

            # per-comparison fold changes
            comparisons = []
            for br in blank_rows:
                ba = float(br["blank_area"]) if not _isnan(br.get("blank_area", float("nan"))) else 0.0
                if ba > 0:
                    fc           = sample_area / ba
                    comp_failed  = fc < fold_threshold
                else:
                    fc           = float("nan")   # no blank signal → always kept
                    comp_failed  = False

                comparisons.append({
                    "blank_name":        br["blank_name"],
                    "blank_area":        ba,
                    "blank_rt":          br.get("blank_rt",  float("nan")),
                    "blank_mz":          br.get("blank_mz",  float("nan")),
                    "mz_delta":          br.get("mz_delta",  float("nan")),
                    "mz_gate_rejected":  bool(br.get("mz_gate_rejected", False)),
                    "fold_change":       fc,
                    "comparison_failed": comp_failed,
                })

            # final decision for this (feature, sample_unit)
            any_failed = any(c["comparison_failed"] for c in comparisons)
            has_blank  = any(not _isnan(c["fold_change"]) for c in comparisons)

            if not has_blank:
                final_decision = "kept"
            elif any_failed:
                if su_type == "all_mean":
                    final_decision = "removed"
                    fids_to_remove.add(fid)
                else:
                    final_decision = "zeroed"
                    for sample in affected_samples:
                        cells_to_zero.append((fid, sample))
            else:
                final_decision = "kept"

            for c in comparisons:
                audit_rows.append({
                    "feature_id":        fid,
                    "sample_unit":       su_name,
                    "unit_type":         su_type,
                    "group":             su_group,
                    "sample_area":       sample_area,
                    "blank_name":        c["blank_name"],
                    "blank_area":        c["blank_area"],
                    "blank_rt":          c["blank_rt"],
                    "blank_mz":          c["blank_mz"],
                    "mz_delta":          c["mz_delta"],
                    "mz_gate_rejected":  c["mz_gate_rejected"],
                    "fold_change":       c["fold_change"],
                    "comparison_failed": c["comparison_failed"],
                    "fold_threshold":    fold_threshold,
                    "decision":          final_decision,
                })

    # --- apply corrections to matrix ---
    mat = matrix.copy()

    if sample_mode == "mean":
        mat = mat.drop(index=list(fids_to_remove), errors="ignore")
    else:
        # zero individual cells
        for fid, sample in cells_to_zero:
            if fid in mat.index and sample in mat.columns:
                mat.loc[fid, sample] = 0.0

        # drop rows where every sample is now 0
        all_zero = (mat == 0).all(axis=1)
        newly_removed = list(mat.index[all_zero])
        mat = mat.loc[~all_zero]
        fids_to_remove.update(newly_removed)

        # upgrade "zeroed" to "removed" for features that ended up fully dropped
        if newly_removed and audit_rows:
            nr_set = set(newly_removed)
            for row in audit_rows:
                if row["feature_id"] in nr_set and row["decision"] == "zeroed":
                    row["decision"] = "removed"

    audit_df = pd.DataFrame(audit_rows)
    return mat, audit_df, list(fids_to_remove)


def _isnan(v):
    try:
        return math.isnan(v)
    except (TypeError, ValueError):
        return v is None


# ---------------------------------------------------------------------------
# Backward-compat: write features_removed_blank.csv
# ---------------------------------------------------------------------------

def _write_removed_log(removed_fids, original_matrix, blank_per_file_df,
                       meta, out_path):
    """
    Write features_removed_blank.csv in the same column format as the old
    pipeline for backward compatibility with any downstream tools.

    Columns: feature_id, mean_rt, mean_mz, mean_sample_area,
             max_blank_area, fold_change
    """
    if not removed_fids:
        pd.DataFrame(columns=["feature_id", "mean_rt", "mean_mz",
                               "mean_sample_area", "max_blank_area",
                               "fold_change"]).to_csv(out_path, index=False)
        return

    removed_set = set(removed_fids)
    rows = []
    max_blank = (blank_per_file_df[blank_per_file_df["feature_id"].isin(removed_set)]
                 .groupby("feature_id")["blank_area"].max())

    for fid in removed_fids:
        if fid not in original_matrix.index:
            continue
        mean_area = float(original_matrix.loc[fid].mean())
        mb        = float(max_blank.get(fid, 0.0))
        fc        = (mean_area / mb) if mb > 0 else float("nan")
        rows.append({
            "feature_id":      fid,
            "mean_sample_area": mean_area,
            "max_blank_area":  mb,
            "fold_change":     fc,
        })

    df = pd.DataFrame(rows).set_index("feature_id")
    if meta is not None and "mean_rt" in meta.columns and "mean_mz" in meta.columns:
        df = df.join(meta[["mean_rt", "mean_mz"]], how="left")
        cols = ["mean_rt", "mean_mz", "mean_sample_area", "max_blank_area", "fold_change"]
        df = df[[c for c in cols if c in df.columns]]

    df.to_csv(out_path)


# ---------------------------------------------------------------------------
# Fallback: reconstruct blank_per_feature from blank_features.csv
# ---------------------------------------------------------------------------

def _fallback_blank_per_file(blank_df):
    """
    Build a synthetic blank_per_feature DataFrame from the legacy
    blank_features.csv (single max-area row per feature).
    Used when blank_per_feature.csv has not yet been generated.
    """
    rows = []
    for fid, row in blank_df.iterrows():
        rows.append({
            "feature_id": fid,
            "blank_name": "BLANK_LEGACY",
            "blank_area": float(row.get("max_blank_area", 0.0)),
            "blank_rt":   row.get("blank_rt",  float("nan")),
            "blank_mz":   row.get("blank_mz",  float("nan")),
        })
    return pd.DataFrame(rows, columns=["feature_id", "blank_name",
                                        "blank_area", "blank_rt", "blank_mz"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 2: blank correction --------------------------------------")

    # --- config parameters (with defaults for backward compat) ---
    sample_mode    = getattr(cfg, "BLANK_SAMPLE_MODE",    "mean")
    ref_mode       = getattr(cfg, "BLANK_REFERENCE_MODE", "max")
    use_mz         = getattr(cfg, "BLANK_USE_MZ",         False)
    mz_tol         = getattr(cfg, "BLANK_MZ_TOLERANCE",   0.005)
    fold_threshold = cfg.FOLD_CHANGE_THRESHOLD

    print(f"  sample mode   : {sample_mode}")
    print(f"  blank ref mode: {ref_mode}")
    print(f"  fold threshold: {fold_threshold}x")
    print(f"  m/z gate      : {'enabled  (BLANK_MZ_TOLERANCE = ' + str(mz_tol) + ' Da)' if use_mz else 'disabled'}")

    # --- load inputs ---
    matrix_path       = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
    blank_pf_path     = os.path.join(cfg.OUTPUT_DIR, "blank_per_feature.csv")
    blank_path        = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")
    metadata_path     = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")
    groups_path       = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")

    for p in (matrix_path, blank_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found - run data_import.py first.")

    matrix   = pd.read_csv(matrix_path, index_col="feature_id")
    blank_df = pd.read_csv(blank_path,  index_col="feature_id")

    print(f"  input   : {matrix.shape[0]} features x {matrix.shape[1]} samples")

    # --- blank per-file table ---
    if os.path.exists(blank_pf_path):
        blank_per_file = pd.read_csv(blank_pf_path)
    else:
        print("  [warning] blank_per_feature.csv not found - re-run data_import.py "
              "for full audit capability.  Falling back to blank_features.csv.")
        blank_per_file = _fallback_blank_per_file(blank_df)

    n_blank_files = blank_per_file["blank_name"].nunique()
    print(f"  blank files   : {n_blank_files} ({', '.join(sorted(blank_per_file['blank_name'].unique()))})")

    # --- feature metadata (needed for m/z gate and audit join) ---
    meta = None
    feat_mz_map = {}
    if os.path.exists(metadata_path):
        meta = pd.read_csv(metadata_path, index_col="feature_id")
        feat_mz_map = meta["mean_mz"].to_dict()

    # --- group map (needed for per_group and per_sample audit labels) ---
    group_map = {}
    if os.path.exists(groups_path):
        gdf       = pd.read_csv(groups_path, index_col="sample")
        group_map = gdf["group"].to_dict()

    # --- apply m/z gate to per-file blank table ---
    if use_mz and feat_mz_map:
        blank_per_file = _apply_mz_gate(blank_per_file, feat_mz_map, mz_tol)
        n_rejected = int(blank_per_file.get("mz_gate_rejected", pd.Series(dtype=bool)).sum())
        if "mz_gate_rejected" in blank_per_file.columns:
            n_rejected = int(blank_per_file["mz_gate_rejected"].sum())
            print(f"  m/z gate      : {n_rejected} (feature, blank) pair(s) rejected "
                  f"(area zeroed; detail in audit log)")
    else:
        if "mz_gate_rejected" not in blank_per_file.columns:
            blank_per_file["mz_gate_rejected"] = False
        if "mz_delta" not in blank_per_file.columns:
            blank_per_file["mz_delta"] = float("nan")

    # --- aggregate blank reference ---
    blank_ref = _aggregate_blank_reference(blank_per_file, ref_mode)

    # --- sample units ---
    sample_units = _compute_sample_units(matrix, group_map, sample_mode)

    # --- build audit and apply corrections ---
    corrected, audit_df, removed_fids = _build_and_apply(
        matrix, blank_ref, sample_units, fold_threshold, sample_mode
    )

    # --- join mean_rt / mean_mz from metadata into audit ---
    if not audit_df.empty and meta is not None:
        audit_df = audit_df.join(
            meta[["mean_rt", "mean_mz"]], on="feature_id", how="left"
        )
        # reorder columns for readability
        front = ["feature_id", "mean_rt", "mean_mz",
                 "sample_unit", "unit_type", "group", "sample_area"]
        back  = [c for c in audit_df.columns if c not in front]
        audit_df = audit_df[front + back]

    # --- keyword filter (name / molecular formula) ---
    exclude_keywords = [kw.lower() for kw in getattr(cfg, "BLANK_EXCLUDE_KEYWORDS", [])]
    if exclude_keywords:
        name_map    = meta["compound_name"].to_dict() if meta is not None and "compound_name" in meta.columns else {}
        formula_map = {}
        enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
        if os.path.exists(enriched_path):
            enriched = pd.read_csv(enriched_path, index_col="feature_id")
            if "molecular_formula" in enriched.columns:
                formula_map = enriched["molecular_formula"].dropna().to_dict()

        kw_removed = []
        for fid in list(corrected.index):
            name    = str(name_map.get(fid, "")).lower()
            formula = str(formula_map.get(fid, "")).lower()
            if any(kw in name or kw in formula for kw in exclude_keywords):
                kw_removed.append(fid)

        if kw_removed:
            corrected = corrected.drop(index=kw_removed, errors="ignore")
            removed_fids = list(set(removed_fids) | set(kw_removed))
            print(f"  keyword filter: {len(kw_removed)} feature(s) removed "
                  f"(keywords: {', '.join(cfg.BLANK_EXCLUDE_KEYWORDS)})")

    # --- count outcomes ---
    n_removed = len(removed_fids)
    n_zeroed_cells = 0
    if not audit_df.empty and "decision" in audit_df.columns:
        # unique (feature, sample_unit) pairs that were zeroed (but feature not removed)
        zeroed_rows = audit_df[audit_df["decision"] == "zeroed"]
        if not zeroed_rows.empty:
            n_zeroed_cells = zeroed_rows.drop_duplicates(
                ["feature_id", "sample_unit"]
            ).shape[0]

    print(f"  removed : {n_removed} feature(s) fully removed")
    if n_zeroed_cells:
        print(f"  zeroed  : {n_zeroed_cells} (feature, sample/group) cell(s) zeroed "
              f"(feature retained for other samples/groups)")
    print(f"  retained: {len(corrected)} features")

    # --- write outputs ---
    out_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    corrected.to_csv(out_path)
    print(f"  -> {out_path}")

    audit_path = os.path.join(cfg.OUTPUT_DIR, "blank_correction_audit.csv")
    if not audit_df.empty:
        audit_df.to_csv(audit_path, index=False)
    else:
        # write empty file with correct headers
        pd.DataFrame(columns=[
            "feature_id", "mean_rt", "mean_mz",
            "sample_unit", "unit_type", "group", "sample_area",
            "blank_name", "blank_area", "blank_rt", "blank_mz",
            "mz_delta", "mz_gate_rejected",
            "fold_change", "comparison_failed", "fold_threshold", "decision",
        ]).to_csv(audit_path, index=False)
    print(f"  -> {audit_path}  ({len(audit_df)} comparison row(s))")

    removed_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_blank.csv")
    _write_removed_log(removed_fids, matrix, blank_per_file, meta, removed_log)
    if removed_fids:
        print(f"  -> {removed_log}  (backward-compat summary, {len(removed_fids)} feature(s))")

    return corrected


if __name__ == "__main__":
    run()
