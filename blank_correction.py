"""
blank_correction.py - Step 2 of the GCMS processing pipeline.

Removes features whose mean sample area does not exceed the maximum blank
area by at least FOLD_CHANGE_THRESHOLD.  Features absent from blanks are
always retained.

Optional m/z gate (BLANK_USE_MZ in config.py)
---------------------------------------------
blank_features.csv stores the m/z of the best RT-matched blank peak alongside
its area.  When BLANK_USE_MZ = True, a blank peak is only accepted as a match
if its m/z is also within BLANK_MZ_TOLERANCE Da of the feature's mean_mz.
Features whose blank peaks fail the m/z check are treated as "not in blank"
and are therefore never removed by the fold-change filter.

Scientific rationale: without m/z gating, a compound that elutes at the same
retention time in both the sample and the blank (but is a chemically different
species with a different mass) would incorrectly trigger removal of a genuine
sample feature.

Note: the EXCLUSION_LIST (biologically relevant features to withhold from
PCA) is applied later in normalization.py, after blank correction, so that
HCA and the volcano plot always operate on the full feature set.

Input  : output/peak_matrix_raw.csv
         output/blank_features.csv      (max_blank_area, blank_rt, blank_mz)
         output/feature_metadata.csv    (mean_rt, mean_mz for audit joins
                                         and the m/z gate when BLANK_USE_MZ=True)
Output : output/peak_matrix_blank_corrected.csv
         output/features_removed_blank.csv       (audit log with fold-change detail)
         output/features_rescued_mz_gate.csv     (only when BLANK_USE_MZ=True and
                                                   any features were rescued)

Usage:
    python blank_correction.py
"""

import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import pandas as pd

import config


# --- Core logic ---------------------------------------------------------------

def blank_correction(matrix, blank_series, fold_change_threshold):
    """
    Filter features from *matrix* based on blank signal.

    Parameters
    ----------
    matrix                : DataFrame  features x samples
    blank_series          : Series     feature_id -> max blank area
                            (already m/z-gated if BLANK_USE_MZ=True, i.e. blank
                            area set to 0 for m/z-mismatched peaks before calling)
    fold_change_threshold : float

    Returns
    -------
    filtered   : DataFrame  features that passed the threshold
    removed_df : DataFrame  audit table for removed features with columns:
                    feature_id, mean_sample_area, max_blank_area, fold_change
    """
    # align blank signal to matrix index; features absent from blanks -> 0
    blank_aligned = blank_series.reindex(matrix.index).fillna(0)
    in_blank      = blank_aligned > 0
    mean_area     = matrix.mean(axis=1)

    # compute fold change only where blank signal exists (avoids divide-by-zero)
    fold_change = mean_area.where(~in_blank).fillna(
        mean_area / blank_aligned.where(blank_aligned > 0)
    )

    # retain feature if it is absent from blanks OR exceeds the fold-change threshold
    keep     = (~in_blank) | (fold_change >= fold_change_threshold)
    filtered = matrix.loc[keep]
    removed  = matrix.index[~keep]

    removed_df = pd.DataFrame({
        "feature_id":       removed,
        "mean_sample_area": mean_area.loc[removed].values,
        "max_blank_area":   blank_aligned.loc[removed].values,
        "fold_change":      fold_change.loc[removed].values,
    })

    return filtered, removed_df


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 2: blank correction --------------------------------------")

    matrix_path   = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
    blank_path    = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")
    metadata_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")

    for p in (matrix_path, blank_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found - run data_import.py first.")

    matrix   = pd.read_csv(matrix_path, index_col="feature_id")
    blank_df = pd.read_csv(blank_path,  index_col="feature_id")

    # blank_series starts as the RT-only matched blank areas
    blank_series = blank_df["max_blank_area"].copy()

    print(f"  input   : {matrix.shape[0]} features x {matrix.shape[1]} samples")
    print(f"  fold-change threshold : {cfg.FOLD_CHANGE_THRESHOLD}x  "
          f"({(blank_series > 0).sum()} features with RT-matched blank signal)")

    # --- optional m/z gate ---------------------------------------------------
    use_mz = getattr(cfg, "BLANK_USE_MZ", False)
    mz_tol = getattr(cfg, "BLANK_MZ_TOLERANCE", 0.005)

    if use_mz:
        print(f"  m/z gate              : enabled  (BLANK_MZ_TOLERANCE = {mz_tol} Da)")
        if "blank_mz" not in blank_df.columns:
            print("  [warning] blank_features.csv has no 'blank_mz' column - "
                  "re-run data_import.py to generate it.  m/z gate skipped.")
        elif not os.path.exists(metadata_path):
            print("  [warning] feature_metadata.csv not found - m/z gate skipped.")
        else:
            meta     = pd.read_csv(metadata_path, index_col="feature_id")
            feat_mz  = meta["mean_mz"].reindex(blank_df.index)
            blank_mz = blank_df["blank_mz"].reindex(blank_df.index)

            has_blank   = blank_series > 0
            mz_ok       = (blank_mz - feat_mz).abs() <= mz_tol
            # mismatch = has an RT-matched blank peak but m/z is too far away
            mz_mismatch = has_blank & ~mz_ok.fillna(False)
            n_rescued   = int(mz_mismatch.sum())

            # zero out blank area for mismatched peaks -> treated as "not in blank"
            blank_series = blank_series.where(~mz_mismatch, 0.0)

            print(f"  m/z gate rescued      : {n_rescued} feature(s) "
                  f"(RT-matched blank peak rejected; |Δm/z| > {mz_tol} Da)")

            if n_rescued > 0:
                rescued_ids = mz_mismatch[mz_mismatch].index
                rescued_df  = pd.DataFrame({
                    "feature_id":         rescued_ids,
                    "mean_mz":            feat_mz.loc[rescued_ids].values,
                    "blank_mz":           blank_mz.loc[rescued_ids].values,
                    "mz_delta":           (blank_mz - feat_mz).abs().loc[rescued_ids].values,
                    "blank_rt":           blank_df["blank_rt"].reindex(rescued_ids).values,
                    "original_blank_area": blank_df["max_blank_area"].reindex(rescued_ids).values,
                })
                # join mean_rt from metadata for context
                if "mean_rt" in meta.columns:
                    rescued_df = (rescued_df.set_index("feature_id")
                                  .join(meta[["mean_rt"]], how="left")
                                  .reset_index())
                    cols = ["feature_id", "mean_rt", "mean_mz",
                            "blank_rt", "blank_mz", "mz_delta", "original_blank_area"]
                    rescued_df = rescued_df[[c for c in cols if c in rescued_df.columns]]

                out_rescued = os.path.join(cfg.OUTPUT_DIR, "features_rescued_mz_gate.csv")
                rescued_df.to_csv(out_rescued, index=False)
                print(f"  -> {out_rescued}  (features retained due to blank m/z mismatch)")
    else:
        print(f"  m/z gate              : disabled  (BLANK_USE_MZ = False)")

    # --- fold-change filter ---------------------------------------------------
    filtered, removed_df = blank_correction(
        matrix, blank_series, cfg.FOLD_CHANGE_THRESHOLD
    )
    print(f"  removed : {len(removed_df)} features  (fold-change filter)")
    print(f"  retained: {len(filtered)} features")

    if len(removed_df) > 0:
        # join mean_rt and mean_mz from feature metadata for full context
        if os.path.exists(metadata_path):
            meta = pd.read_csv(metadata_path, index_col="feature_id")[["mean_rt", "mean_mz"]]
            removed_df = (removed_df.set_index("feature_id")
                          .join(meta, how="left")
                          .reset_index())
            cols = ["feature_id", "mean_rt", "mean_mz",
                    "mean_sample_area", "max_blank_area", "fold_change"]
            removed_df = removed_df[[c for c in cols if c in removed_df.columns]]

        blank_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_blank.csv")
        removed_df.to_csv(blank_log, index=False)
        print(f"  -> {blank_log}  (incl. fold-change detail)")

    out_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    filtered.to_csv(out_path)
    print(f"  -> {out_path}")

    return filtered


if __name__ == "__main__":
    run()
