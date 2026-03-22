"""
blank_correction.py - Step 2 of the GCMS processing pipeline.

Applies two sequential filters to the raw peak matrix:

  1. Blank correction - removes features whose mean sample area does not exceed
     the maximum blank area by at least FOLD_CHANGE_THRESHOLD.
     Features absent from blanks are always retained.

  2. Exclusion list - removes features whose mean RT falls within
     +-EXCLUSION_RT_MARGIN of any retention time listed in EXCLUSION_LIST.
     Use this to strip known interference compounds (solvent artefacts,
     column bleed, etc.) regardless of their blank signal.

Both filters are applied in one pass; a combined audit log is written.

Input  : output/peak_matrix_raw.csv
         output/blank_features.csv
         output/feature_metadata.csv
Output : output/peak_matrix_blank_corrected.csv
         output/features_removed_blank.csv     (blank filter audit log)
         output/features_removed_exclusion.csv (exclusion list audit log)

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

def apply_exclusion_list(matrix, metadata, exclusion_rts, rt_margin):
    """
    Remove features whose mean RT is within +-rt_margin of any RT in
    *exclusion_rts*.

    Parameters
    ----------
    matrix        : DataFrame  features x samples
    metadata      : DataFrame  feature_id -> mean_rt, mean_mz
    exclusion_rts : list of float  retention times to exclude (minutes)
    rt_margin     : float          +- window in minutes

    Returns
    -------
    filtered : DataFrame  features that passed the exclusion filter
    removed  : DataFrame  excluded features with columns:
                          feature_id, mean_rt, matched_exclusion_rt
    """
    if not exclusion_rts:
        return matrix, pd.DataFrame(columns=["feature_id", "mean_rt",
                                              "matched_exclusion_rt"])

    rt_lookup  = metadata["mean_rt"].reindex(matrix.index)
    hit_rt     = {}   # feature_id -> the exclusion RT that triggered removal

    for fid, feat_rt in rt_lookup.items():
        if pd.isna(feat_rt):
            continue
        for excl_rt in exclusion_rts:
            if abs(feat_rt - excl_rt) <= rt_margin:
                hit_rt[fid] = excl_rt
                break          # one match is enough to exclude

    exclude_idx = list(hit_rt.keys())
    filtered    = matrix.drop(index=exclude_idx)
    removed     = pd.DataFrame([
        {"feature_id": fid, "mean_rt": rt_lookup[fid],
         "matched_exclusion_rt": hit_rt[fid]}
        for fid in exclude_idx
    ])

    return filtered, removed


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

    print("-- Step 2: blank correction & exclusion list ---------------------")

    matrix_path   = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
    blank_path    = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")
    metadata_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")

    for p in (matrix_path, blank_path, metadata_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found - run data_import.py first.")

    matrix       = pd.read_csv(matrix_path,   index_col="feature_id")
    blank_series = pd.read_csv(blank_path,     index_col="feature_id")["max_blank_area"]
    metadata     = pd.read_csv(metadata_path,  index_col="feature_id")

    print(f"  input   : {matrix.shape[0]} features x {matrix.shape[1]} samples")

    # -- 1. blank correction ---------------------------------------------------
    print(f"  blank correction  - fold-change threshold: {cfg.FOLD_CHANGE_THRESHOLD}x  "
          f"({(blank_series > 0).sum()} features with blank signal)")

    filtered, removed_blank = blank_correction(
        matrix, blank_series, cfg.FOLD_CHANGE_THRESHOLD
    )
    print(f"    removed {len(removed_blank)} features")

    if len(removed_blank) > 0:
        blank_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_blank.csv")
        pd.DataFrame({"feature_id": removed_blank}).to_csv(blank_log, index=False)
        print(f"    -> {blank_log}")

    # -- 2. exclusion list -----------------------------------------------------
    excl_rts = [rt for rt in cfg.EXCLUSION_LIST if rt is not None]
    print(f"  exclusion list    - {len(excl_rts)} RT(s) listed, "
          f"margin +-{cfg.EXCLUSION_RT_MARGIN} min")

    filtered, removed_excl = apply_exclusion_list(
        filtered, metadata, excl_rts, cfg.EXCLUSION_RT_MARGIN
    )
    print(f"    removed {len(removed_excl)} features")

    if len(removed_excl) > 0:
        excl_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_exclusion.csv")
        removed_excl.to_csv(excl_log, index=False)
        print(f"    -> {excl_log}")

    # -- summary & output ------------------------------------------------------
    total_removed = len(removed_blank) + len(removed_excl)
    print(f"  retained: {len(filtered)} features  "
          f"(removed {total_removed} total)")

    out_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    filtered.to_csv(out_path)
    print(f"  -> {out_path}")

    return filtered


if __name__ == "__main__":
    run()
