"""
blank_correction.py - Step 2 of the GCMS processing pipeline.

Removes features whose mean sample area does not exceed the maximum blank
area by at least FOLD_CHANGE_THRESHOLD.  Features absent from blanks are
always retained.

Note: the EXCLUSION_LIST (biologically relevant features to withhold from
PCA) is applied later in normalization.py, after blank correction, so that
HCA and the volcano plot always operate on the full feature set.

Input  : output/peak_matrix_raw.csv
         output/blank_features.csv
Output : output/peak_matrix_blank_corrected.csv
         output/features_removed_blank.csv     (audit log)

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
    fold_change_threshold : float

    Returns
    -------
    filtered : DataFrame  features that passed the threshold
    removed  : Index      features that were removed
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

    return filtered, removed


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 2: blank correction --------------------------------------")

    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
    blank_path  = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")

    for p in (matrix_path, blank_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found - run data_import.py first.")

    matrix       = pd.read_csv(matrix_path, index_col="feature_id")
    blank_series = pd.read_csv(blank_path,  index_col="feature_id")["max_blank_area"]

    print(f"  input   : {matrix.shape[0]} features x {matrix.shape[1]} samples")
    print(f"  fold-change threshold: {cfg.FOLD_CHANGE_THRESHOLD}x  "
          f"({(blank_series > 0).sum()} features with blank signal)")

    filtered, removed_blank = blank_correction(
        matrix, blank_series, cfg.FOLD_CHANGE_THRESHOLD
    )
    print(f"  removed : {len(removed_blank)} features")
    print(f"  retained: {len(filtered)} features")

    if len(removed_blank) > 0:
        blank_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_blank.csv")
        pd.DataFrame({"feature_id": removed_blank}).to_csv(blank_log, index=False)
        print(f"  -> {blank_log}")

    out_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    filtered.to_csv(out_path)
    print(f"  -> {out_path}")

    return filtered


if __name__ == "__main__":
    run()
