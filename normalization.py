"""
normalization.py - Step 3 of the GCMS processing pipeline.

Applies three sequential transformations and produces two output matrices:

  1. Normalization   - corrects for differences in total sample signal
                       caused by variation in injection volume or extraction yield.
                       Method: sum normalization (each sample divided by its total
                       signal, rescaled to the median column sum so absolute
                       magnitudes are preserved).

  2. Log2 transformation - compresses the dynamic range, stabilises variance
                       across features, and brings the distribution closer to
                       normal. log2(x + 1) handles zero values gracefully.

  3. Pareto scaling  - mean-centres each feature and divides by its square-root
                       standard deviation. This reduces the dominance of
                       high-intensity features while preserving more biological
                       variance structure than unit-variance (auto) scaling.
                       Recommended for untargeted GC-MS metabolomics.

Zero-variance features (identical in all samples) are automatically dropped
before scaling, as they carry no discriminating information.

Output orientation: samples x features  (rows = samples, columns = features)
This is the standard input format for scikit-learn PCA, HCA, OPLS-DA, etc.

Two output matrices are produced:

  peak_matrix_processed.csv       full feature set - used by HCA and volcano plot.
                                  No exclusion list applied so that all biologically
                                  relevant features remain visible in those analyses.

  peak_matrix_processed_pca.csv   exclusion-list filtered - used by PCA only.
                                  Features listed in EXCLUSION_LIST are removed
                                  before normalization so that already-characterised
                                  compounds do not dominate the principal components.
                                  Identical to peak_matrix_processed.csv when
                                  EXCLUSION_LIST is empty.

Input  : output/peak_matrix_blank_corrected.csv
         output/feature_metadata.csv   (for exclusion list RT lookup)
Output : output/peak_matrix_processed.csv
         output/peak_matrix_processed_pca.csv
         output/features_removed_exclusion.csv  (audit log, only if list is non-empty)

Usage:
    python normalization.py
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

import config


# --- Prevalence filter -------------------------------------------------------

def _prevalence_filter(matrix, min_prevalence):
    """
    Remove features detected in fewer than *min_prevalence* fraction of samples.
    'Detected' means area > 0 (operates on a features x samples matrix).

    Parameters
    ----------
    matrix          : DataFrame  features x samples  (values are raw areas)
    min_prevalence  : float      0.0 = keep all; 1.0 = require detection in every sample

    Returns
    -------
    filtered : DataFrame
    n_removed : int
    """
    if min_prevalence <= 0.0:
        return matrix, 0
    detected_fraction = (matrix > 0).mean(axis=1)
    keep = detected_fraction >= min_prevalence
    return matrix.loc[keep], int((~keep).sum())


# --- Exclusion list (PCA only) -----------------------------------------------

def _apply_exclusion(matrix, metadata, exclusion_rts, rt_margin):
    """
    Remove features whose mean RT falls within +-rt_margin of any RT in
    *exclusion_rts*.  Returns (filtered_matrix, removed_df).
    Called only when building the PCA-specific matrix.
    """
    if not exclusion_rts:
        return matrix, None

    rt_lookup = metadata["mean_rt"].reindex(matrix.index)
    hit_rt    = {}

    for fid, feat_rt in rt_lookup.items():
        if pd.isna(feat_rt):
            continue
        for excl_rt in exclusion_rts:
            if abs(feat_rt - excl_rt) <= rt_margin:
                hit_rt[fid] = excl_rt
                break

    exclude_idx = list(hit_rt.keys())
    filtered    = matrix.drop(index=exclude_idx)
    removed     = pd.DataFrame([
        {"feature_id": fid, "mean_rt": float(rt_lookup[fid]),
         "matched_exclusion_rt": hit_rt[fid]}
        for fid in exclude_idx
    ])
    return filtered, removed


# --- Normalization methods ----------------------------------------------------

def sum_normalize(matrix):
    """
    Divide each sample (column) by its total signal, then rescale by the
    median column sum so that values remain in a biologically interpretable range.

    Operating on features x samples matrix (column-wise normalization).
    """
    col_sums      = matrix.sum(axis=0)
    median_sum    = col_sums.median()
    if median_sum == 0:
        raise ValueError("All sample sums are zero - check input data.")
    return matrix.divide(col_sums, axis=1) * median_sum


def median_normalize(matrix):
    """
    Scale each sample so its median feature value equals the global median.
    More robust than sum normalization when many features are zero.
    """
    col_medians   = matrix.median(axis=0)
    global_median = col_medians.median()
    if global_median == 0:
        raise ValueError("Global median is zero - check input data.")
    return matrix.divide(col_medians.replace(0, global_median), axis=1) * global_median


# --- Log transformation -------------------------------------------------------

def log_transform(matrix, base=2):
    """
    Apply log(x + 1) transformation to compress dynamic range.

    Parameters
    ----------
    base : int or float  - 2 (log2) or math.e (natural log)
    """
    if base == 2:
        return np.log2(matrix + 1)
    elif base == math.e or base == "e":
        return np.log1p(matrix)
    else:
        return np.log(matrix + 1) / np.log(base)


# --- Scaling methods ---------------------------------------------------------

def pareto_scale(matrix):
    """
    Pareto scaling: (x - mean_feature) / sqrtstd_feature

    Applied row-wise (per feature) across all samples.
    Features with zero variance after log transform are dropped.
    """
    std_vec  = matrix.std(axis=1, ddof=1)
    zero_var = std_vec == 0
    if zero_var.any():
        n = zero_var.sum()
        print(f"  [info] dropping {n} zero-variance feature(s) before scaling")
        matrix  = matrix.loc[~zero_var]
        std_vec = std_vec.loc[~zero_var]

    mean_vec  = matrix.mean(axis=1)
    sqrt_std  = std_vec.apply(math.sqrt)
    return matrix.subtract(mean_vec, axis=0).divide(sqrt_std, axis=0)


def auto_scale(matrix):
    """
    Auto (unit-variance) scaling: (x - mean_feature) / std_feature

    Gives every feature equal variance - stronger equalisation than Pareto,
    but may amplify noise in low-abundance features.
    """
    std_vec  = matrix.std(axis=1, ddof=1)
    zero_var = std_vec == 0
    if zero_var.any():
        n = zero_var.sum()
        print(f"  [info] dropping {n} zero-variance feature(s) before scaling")
        matrix  = matrix.loc[~zero_var]
        std_vec = std_vec.loc[~zero_var]

    mean_vec = matrix.mean(axis=1)
    return matrix.subtract(mean_vec, axis=0).divide(std_vec, axis=0)


# --- Shared transformation helper --------------------------------------------

def _transform(matrix, cfg):
    """Apply normalization + log transform + scaling to a features x samples matrix."""
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)

    if cfg.NORMALIZATION == "sum":
        matrix = sum_normalize(matrix)
    elif cfg.NORMALIZATION == "median":
        matrix = median_normalize(matrix)

    matrix = log_transform(matrix, cfg.LOG_BASE)

    if cfg.SCALING == "pareto":
        matrix = pareto_scale(matrix)
    elif cfg.SCALING == "auto":
        matrix = auto_scale(matrix)

    processed              = matrix.T
    processed.index.name   = "sample"
    processed.columns.name = "feature_id"
    return processed


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 3: normalization -----------------------------------------")

    in_path       = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    metadata_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")

    for p in (in_path, metadata_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found - run data_import.py first.")

    matrix   = pd.read_csv(in_path,       index_col="feature_id")
    metadata = pd.read_csv(metadata_path, index_col="feature_id")

    base_label = "e" if cfg.LOG_BASE == math.e else cfg.LOG_BASE
    print(f"  input  : {matrix.shape[0]} features x {matrix.shape[1]} samples")
    print(f"  steps  : {cfg.NORMALIZATION} normalization  ->  "
          f"log{base_label}(x+1)  ->  {cfg.SCALING} scaling")

    # -- full matrix (HCA + volcano) ------------------------------------------
    matrix_hca, n_rem_hca = _prevalence_filter(matrix, cfg.MIN_PREVALENCE_HCA)
    if n_rem_hca:
        print(f"  prevalence filter (HCA): removed {n_rem_hca} feature(s) "
              f"detected in < {cfg.MIN_PREVALENCE_HCA*100:.0f}% of samples")

    processed = _transform(matrix_hca, cfg)
    out_full  = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed.csv")
    processed.to_csv(out_full)
    print(f"  -> {out_full}  ({processed.shape[0]} samples x {processed.shape[1]} features)"
          f"  [HCA, volcano]")

    # -- exclusion-filtered + prevalence-filtered matrix (PCA only) -----------
    matrix_pca, n_rem_pca = _prevalence_filter(matrix, cfg.MIN_PREVALENCE_PCA)
    if n_rem_pca:
        print(f"  prevalence filter (PCA): removed {n_rem_pca} feature(s) "
              f"detected in < {cfg.MIN_PREVALENCE_PCA*100:.0f}% of samples")

    excl_rts = [rt for rt in cfg.EXCLUSION_LIST if rt is not None]

    if excl_rts:
        matrix_pca, removed = _apply_exclusion(
            matrix_pca, metadata, excl_rts, cfg.EXCLUSION_RT_MARGIN
        )
        print(f"  exclusion list : {len(excl_rts)} RT(s), "
              f"margin +-{cfg.EXCLUSION_RT_MARGIN} min  ->  "
              f"removed {len(removed)} feature(s) for PCA")

        excl_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_exclusion.csv")
        removed.to_csv(excl_log, index=False)
        print(f"  -> {excl_log}")

    processed_pca = _transform(matrix_pca, cfg)
    if not excl_rts and n_rem_pca == 0:
        processed_pca = processed   # identical when no filters applied

    out_pca = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_pca.csv")
    processed_pca.to_csv(out_pca)
    print(f"  -> {out_pca}  ({processed_pca.shape[0]} samples x {processed_pca.shape[1]} features)"
          f"  [PCA]")

    return processed, processed_pca


if __name__ == "__main__":
    run()
