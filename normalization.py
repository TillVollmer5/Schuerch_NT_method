"""
normalization.py - Step 3 of the GCMS processing pipeline.

Applies three sequential transformations to produce a matrix suitable for
principal component analysis (PCA) and other multivariate analyses:

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
This is the standard input format for scikit-learn PCA, OPLS-DA, etc.

Input  : output/peak_matrix_blank_corrected.csv
Output : output/peak_matrix_processed.csv

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


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 3: normalization -----------------------------------------")

    in_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    if not os.path.exists(in_path):
        raise FileNotFoundError(
            f"{in_path} not found - run blank_correction.py first."
        )

    matrix = pd.read_csv(in_path, index_col="feature_id")
    print(f"  input  : {matrix.shape[0]} features x {matrix.shape[1]} samples")

    # convert all columns to numeric (safety check)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)

    # -- 1. normalization ------------------------------------------------------
    if cfg.NORMALIZATION == "sum":
        matrix = sum_normalize(matrix)
        print("  step 1 : sum normalization")
    elif cfg.NORMALIZATION == "median":
        matrix = median_normalize(matrix)
        print("  step 1 : median normalization")
    else:
        print("  step 1 : normalization skipped")

    # -- 2. log transformation -------------------------------------------------
    matrix = log_transform(matrix, cfg.LOG_BASE)
    base_label = "e" if cfg.LOG_BASE == math.e else cfg.LOG_BASE
    print(f"  step 2 : log{base_label}(x + 1) transformation")

    # -- 3. scaling ------------------------------------------------------------
    if cfg.SCALING == "pareto":
        matrix = pareto_scale(matrix)
        print("  step 3 : Pareto scaling")
    elif cfg.SCALING == "auto":
        matrix = auto_scale(matrix)
        print("  step 3 : auto (unit-variance) scaling")
    else:
        print("  step 3 : scaling skipped")

    # -- output: transpose to samples x features for multivariate analysis -----
    processed              = matrix.T
    processed.index.name   = "sample"
    processed.columns.name = "feature_id"

    out_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed.csv")
    processed.to_csv(out_path)
    print(f"  -> {out_path}")
    print(f"     {processed.shape[0]} samples x {processed.shape[1]} features")
    print("     orientation: samples x features  (ready for PCA / OPLS-DA)")

    return processed


if __name__ == "__main__":
    run()
