"""
normalization.py - Step 3 of the GCMS processing pipeline.

Applies three sequential transformations and produces three output matrices,
each with settings optimised for a different downstream analysis:

  peak_matrix_processed_pca.csv   PCA    - exclusion-list filtered, per-analysis settings
  peak_matrix_processed_hca.csv   HCA    - full feature set, per-analysis settings
  peak_matrix_processed_volcano.csv  Volcano - full feature set, no scaling

Per-analysis normalization / log / scaling methods are chosen via the
NORMALIZATION_<X>, LOG_BASE_<X>, and SCALING_<X> keys in config.py, falling
back to the global NORMALIZATION, LOG_BASE, and SCALING when the override is None.

Normalization methods (sample-wise, on raw areas before log):
  "none"   - skip
  "sum"    - divide by total signal, rescale to median column sum
  "median" - scale so sample median equals global median
  "pqn"    - Probabilistic Quotient Normalization (Dieterle et al. 2006)

Log transform methods (applied after normalization):
  2        - log2(x+1)
  10       - log10(x+1)
  math.e   - natural log ln(x+1)
  "sqrt"   - square-root transform sqrt(x)
  "cbrt"   - cube-root transform cbrt(x)
  "none"   - skip

Scaling methods (feature-wise, across all samples, after log):
  "none"   - skip
  "pareto" - (x - mean) / sqrt(std)
  "auto"   - (x - mean) / std  (unit-variance)
  "vast"   - (x - mean) * mean / std²  (VAST; Van den Berg et al. 2006)
  "range"  - (x - mean) / (max - min)
  "level"  - (x - mean) / mean

Output orientation: samples x features  (rows = samples, columns = features)
This is the standard input format for scikit-learn PCA, seaborn clustermap, etc.

Audit log files written (when applicable):
  features_removed_exclusion.csv       - exclusion list matches (PCA only)
  features_removed_prevalence_hca.csv  - sparse features dropped for HCA
  features_removed_prevalence_pca.csv  - sparse features dropped for PCA
  features_removed_prevalence_volcano.csv - sparse features dropped for volcano
  features_removed_zero_variance.csv   - zero-variance features dropped by scaling

Input  : output/peak_matrix_blank_corrected.csv
         output/feature_metadata.csv
Output : output/peak_matrix_processed_pca.csv
         output/peak_matrix_processed_hca.csv
         output/peak_matrix_processed_volcano.csv

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


# --- Analysis settings helper ------------------------------------------------

def _get_settings(cfg, analysis):
    """
    Return (normalization, log_base, scaling) for a named analysis,
    falling back to the global defaults when the override is None.

    analysis : "pca" | "hca" | "volcano" | "boxplot"
    """
    key = analysis.upper()
    norm    = getattr(cfg, f"NORMALIZATION_{key}", None) or cfg.NORMALIZATION
    log_b   = getattr(cfg, f"LOG_BASE_{key}",     None)
    scaling = getattr(cfg, f"SCALING_{key}",       None) or cfg.SCALING

    # LOG_BASE override: None means inherit global; "none" means skip
    if log_b is None:
        log_b = cfg.LOG_BASE

    return norm, log_b, scaling


# --- Prevalence filter -------------------------------------------------------

def _prevalence_filter(matrix, min_prevalence, metadata=None):
    """
    Remove features detected in fewer than *min_prevalence* fraction of samples.
    'Detected' means area > 0 (operates on a features x samples matrix).
    """
    if min_prevalence <= 0.0:
        return matrix, pd.DataFrame()

    detected_fraction = (matrix > 0).mean(axis=1)
    n_detected        = (matrix > 0).sum(axis=1)
    keep              = detected_fraction >= min_prevalence
    removed_ids       = matrix.index[~keep]

    removed_df = pd.DataFrame({
        "feature_id":         removed_ids,
        "n_samples_detected": n_detected.loc[removed_ids].values,
        "n_samples_total":    matrix.shape[1],
        "detected_fraction":  detected_fraction.loc[removed_ids].values,
    })

    if metadata is not None and not removed_df.empty:
        meta_cols = [c for c in ("mean_rt", "mean_mz") if c in metadata.columns]
        if meta_cols:
            removed_df = (removed_df.set_index("feature_id")
                          .join(metadata[meta_cols], how="left")
                          .reset_index())
            cols = ["feature_id"] + meta_cols + ["n_samples_detected",
                                                   "n_samples_total",
                                                   "detected_fraction"]
            removed_df = removed_df[[c for c in cols if c in removed_df.columns]]

    return matrix.loc[keep], removed_df


# --- Exclusion list (PCA only) -----------------------------------------------

def _apply_exclusion(matrix, metadata, exclusion_entries, rt_margin, mz_tolerance):
    """
    Remove features whose mean_rt and mean_mz fall within the given tolerances
    of any entry in *exclusion_entries*.  Returns (filtered_matrix, removed_df).
    """
    if not exclusion_entries:
        return matrix, None

    rt_lookup = metadata["mean_rt"].reindex(matrix.index)
    mz_lookup = metadata["mean_mz"].reindex(matrix.index) if "mean_mz" in metadata.columns else None
    hit_rt    = {}
    hit_dev   = {}

    for fid, feat_rt in rt_lookup.items():
        if pd.isna(feat_rt):
            continue
        feat_mz = float(mz_lookup[fid]) if (mz_lookup is not None and not pd.isna(mz_lookup[fid])) else None
        for (excl_rt, excl_mz) in exclusion_entries:
            if abs(feat_rt - excl_rt) > rt_margin:
                continue
            if excl_mz is not None and feat_mz is not None:
                if abs(feat_mz - excl_mz) > mz_tolerance:
                    continue
            hit_rt[fid]  = excl_rt
            hit_dev[fid] = round(abs(feat_rt - excl_rt), 6)
            break

    exclude_idx = list(hit_rt.keys())
    filtered    = matrix.drop(index=exclude_idx)
    removed     = pd.DataFrame([
        {
            "feature_id":           fid,
            "mean_rt":              float(rt_lookup[fid]),
            "matched_exclusion_rt": hit_rt[fid],
            "rt_deviation":         hit_dev[fid],
        }
        for fid in exclude_idx
    ])
    return filtered, removed


# --- Normalization methods ----------------------------------------------------

def sum_normalize(matrix):
    """
    Divide each sample (column) by its total signal, then rescale by the
    median column sum so that values remain in a biologically interpretable range.
    """
    col_sums   = matrix.sum(axis=0)
    median_sum = col_sums.median()
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


def pqn_normalize(matrix):
    """
    Probabilistic Quotient Normalization (Dieterle et al. 2006, Anal. Chem.).

    Algorithm:
      1. Compute reference spectrum as the feature-wise median across all samples.
      2. For each sample, divide each detected feature value by the reference.
         (Features where reference = 0 are excluded from quotient calculation.)
      3. The PQN coefficient = median of those quotients for the sample.
      4. Divide the sample by its PQN coefficient.
      5. Rescale by the median of all PQN coefficients to preserve scale.

    This corrects for global dilution/concentration differences while preserving
    genuine biological fold-changes — unlike sum normalization which deflates
    all features when total signal differs biologically.
    """
    matrix    = matrix.copy().apply(pd.to_numeric, errors="coerce").fillna(0)
    reference = matrix.median(axis=1)      # feature-wise median = reference spectrum

    pqn_factors = []
    for col in matrix.columns:
        sample    = matrix[col]
        mask      = reference > 0          # only use features with non-zero reference
        quotients = sample[mask] / reference[mask]
        factor    = float(quotients.median())
        if factor <= 0 or np.isnan(factor):
            factor = 1.0
        pqn_factors.append(factor)
        matrix[col] = sample / factor

    # rescale so overall magnitude is preserved (multiply by median factor)
    median_factor = float(np.median(pqn_factors))
    if median_factor > 0:
        matrix = matrix * median_factor

    return matrix


# --- Log transform ------------------------------------------------------------

def log_transform(matrix, base=2):
    """
    Apply log(x + 1) or power-root transform to compress dynamic range.

    base : 2       -> log2(x+1)
           10      -> log10(x+1)
           math.e  -> ln(x+1)
           "e"     -> ln(x+1)
           "sqrt"  -> sqrt(x)
           "cbrt"  -> cbrt(x)
           "none"  -> identity (no transform)
    """
    if base == "none" or base is None:
        return matrix
    if base == "sqrt":
        return np.sqrt(matrix.clip(lower=0))
    if base == "cbrt":
        return matrix.apply(np.cbrt)
    if base == 2:
        return np.log2(matrix + 1)
    if base == 10:
        return np.log10(matrix + 1)
    if base == math.e or base == "e":
        return np.log1p(matrix)
    # generic base
    return np.log(matrix + 1) / np.log(base)


# --- Scaling methods ---------------------------------------------------------

def _drop_bad(matrix, mask, label):
    removed = matrix.index[mask].tolist()
    if mask.any():
        n = mask.sum()
        print(f"  [info] dropping {n} feature(s) with zero {label} before scaling")
        matrix = matrix.loc[~mask]
    return matrix, removed


def pareto_scale(matrix):
    """Pareto: (x - mean) / sqrt(std)"""
    std_vec  = matrix.std(axis=1, ddof=1)
    matrix, removed = _drop_bad(matrix, std_vec == 0, "variance")
    std_vec  = matrix.std(axis=1, ddof=1)
    mean_vec = matrix.mean(axis=1)
    sqrt_std = std_vec.apply(math.sqrt)
    return matrix.subtract(mean_vec, axis=0).divide(sqrt_std, axis=0), removed


def auto_scale(matrix):
    """Auto (unit-variance): (x - mean) / std"""
    std_vec  = matrix.std(axis=1, ddof=1)
    matrix, removed = _drop_bad(matrix, std_vec == 0, "variance")
    std_vec  = matrix.std(axis=1, ddof=1)
    mean_vec = matrix.mean(axis=1)
    return matrix.subtract(mean_vec, axis=0).divide(std_vec, axis=0), removed


def vast_scale(matrix):
    """
    VAST (Variable Stability) scaling: (x - mean) * mean / std²
    Equivalent to auto-scaling weighted by 1/CV.
    Favours stable (low CV) features; de-emphasises noisy ones.
    (Van den Berg et al. 2006, BMC Genomics)
    """
    std_vec  = matrix.std(axis=1, ddof=1)
    mean_vec = matrix.mean(axis=1)
    bad      = (std_vec == 0) | (mean_vec == 0)
    matrix, removed = _drop_bad(matrix, bad, "variance or mean")
    std_vec  = matrix.std(axis=1, ddof=1)
    mean_vec = matrix.mean(axis=1)
    denom    = std_vec ** 2 / mean_vec       # std² / mean
    return matrix.subtract(mean_vec, axis=0).divide(denom, axis=0), removed


def range_scale(matrix):
    """Range scaling: (x - mean) / (max - min)"""
    mean_vec  = matrix.mean(axis=1)
    range_vec = matrix.max(axis=1) - matrix.min(axis=1)
    matrix, removed = _drop_bad(matrix, range_vec == 0, "range")
    mean_vec  = matrix.mean(axis=1)
    range_vec = matrix.max(axis=1) - matrix.min(axis=1)
    return matrix.subtract(mean_vec, axis=0).divide(range_vec, axis=0), removed


def level_scale(matrix):
    """Level scaling: (x - mean) / mean"""
    mean_vec = matrix.mean(axis=1)
    matrix, removed = _drop_bad(matrix, mean_vec == 0, "mean")
    mean_vec = matrix.mean(axis=1)
    return matrix.subtract(mean_vec, axis=0).divide(mean_vec, axis=0), removed


# --- Shared transformation helper --------------------------------------------

def _transform(matrix, normalization, log_base, scaling):
    """
    Apply normalization + log transform + scaling to a features x samples matrix.

    Returns
    -------
    processed    : DataFrame  samples x features (transposed)
    zero_var_ids : list of str  feature_ids dropped by scaling (zero variance)
    """
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)

    # --- normalization (sample-wise, on raw areas) ---------------------------
    if normalization == "sum":
        matrix = sum_normalize(matrix)
    elif normalization == "median":
        matrix = median_normalize(matrix)
    elif normalization == "pqn":
        matrix = pqn_normalize(matrix)
    # "none" or anything else: skip

    # --- log transform -------------------------------------------------------
    matrix = log_transform(matrix, log_base)

    # --- scaling (feature-wise, after log) -----------------------------------
    zero_var_ids = []
    if scaling == "pareto":
        matrix, zero_var_ids = pareto_scale(matrix)
    elif scaling == "auto":
        matrix, zero_var_ids = auto_scale(matrix)
    elif scaling == "vast":
        matrix, zero_var_ids = vast_scale(matrix)
    elif scaling == "range":
        matrix, zero_var_ids = range_scale(matrix)
    elif scaling == "level":
        matrix, zero_var_ids = level_scale(matrix)
    # "none" or anything else: skip

    processed              = matrix.T
    processed.index.name   = "sample"
    processed.columns.name = "feature_id"
    return processed, zero_var_ids


def _label(log_base):
    """Human-readable log label for console output."""
    if log_base == "none" or log_base is None:
        return "none"
    if log_base == math.e or log_base == "e":
        return "ln"
    if log_base == "sqrt":
        return "sqrt"
    if log_base == "cbrt":
        return "cbrt"
    return f"log{log_base}"


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

    print(f"  input  : {matrix.shape[0]} features x {matrix.shape[1]} samples")

    zero_var_log = {}   # analysis -> list of dropped feature_ids

    # =========================================================================
    # HCA matrix  (full feature set, HCA-specific normalization/log/scaling)
    # =========================================================================
    norm_hca, log_hca, scal_hca = _get_settings(cfg, "hca")
    print(f"  HCA    : {norm_hca} norm  ->  {_label(log_hca)}(x+1)  ->  {scal_hca} scaling")

    matrix_hca, removed_prev_hca = _prevalence_filter(
        matrix, cfg.MIN_PREVALENCE_HCA, metadata
    )
    if not removed_prev_hca.empty:
        print(f"  prevalence filter (HCA): removed {len(removed_prev_hca)} feature(s) "
              f"detected in < {cfg.MIN_PREVALENCE_HCA*100:.0f}% of samples")
        out = os.path.join(cfg.OUTPUT_DIR, "features_removed_prevalence_hca.csv")
        removed_prev_hca.to_csv(out, index=False)
        print(f"  -> {out}")

    processed_hca, zv_hca = _transform(matrix_hca, norm_hca, log_hca, scal_hca)
    zero_var_log["hca"] = zv_hca

    out_hca = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_hca.csv")
    processed_hca.to_csv(out_hca)
    print(f"  -> {out_hca}  ({processed_hca.shape[0]} samples x {processed_hca.shape[1]} features)")

    # =========================================================================
    # Volcano matrix (full feature set, no scaling — required for fold-change)
    # =========================================================================
    norm_vol, log_vol, scal_vol = _get_settings(cfg, "volcano")
    if scal_vol not in ("none", None, ""):
        print(f"  [warning] SCALING_VOLCANO = '{scal_vol}' will distort fold-changes. "
              f"Setting to 'none' for volcano.")
        scal_vol = "none"
    print(f"  Volcano: {norm_vol} norm  ->  {_label(log_vol)}(x+1)  ->  none scaling")

    min_prev_vol = getattr(cfg, "MIN_PREVALENCE_VOLCANO", 0.0)
    matrix_vol, removed_prev_vol = _prevalence_filter(
        matrix, min_prev_vol, metadata
    )
    if not removed_prev_vol.empty:
        print(f"  prevalence filter (Volcano): removed {len(removed_prev_vol)} feature(s) "
              f"detected in < {min_prev_vol*100:.0f}% of samples")
        out = os.path.join(cfg.OUTPUT_DIR, "features_removed_prevalence_volcano.csv")
        removed_prev_vol.to_csv(out, index=False)
        print(f"  -> {out}")

    processed_vol, _ = _transform(matrix_vol, norm_vol, log_vol, "none")

    out_vol = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_volcano.csv")
    processed_vol.to_csv(out_vol)
    print(f"  -> {out_vol}  ({processed_vol.shape[0]} samples x {processed_vol.shape[1]} features)")

    # =========================================================================
    # PCA matrix  (exclusion-list filtered, PCA-specific settings)
    # =========================================================================
    norm_pca, log_pca, scal_pca = _get_settings(cfg, "pca")
    print(f"  PCA    : {norm_pca} norm  ->  {_label(log_pca)}(x+1)  ->  {scal_pca} scaling")

    matrix_pca, removed_prev_pca = _prevalence_filter(
        matrix, cfg.MIN_PREVALENCE_PCA, metadata
    )
    if not removed_prev_pca.empty:
        print(f"  prevalence filter (PCA): removed {len(removed_prev_pca)} feature(s) "
              f"detected in < {cfg.MIN_PREVALENCE_PCA*100:.0f}% of samples")
        out = os.path.join(cfg.OUTPUT_DIR, "features_removed_prevalence_pca.csv")
        removed_prev_pca.to_csv(out, index=False)
        print(f"  -> {out}")

    excl_entries = [(e[0], e[1]) for e in cfg.EXCLUSION_LIST
                    if e[0] is not None and e[1] is not None]
    if excl_entries:
        matrix_pca, removed_excl = _apply_exclusion(
            matrix_pca, metadata, excl_entries,
            cfg.EXCLUSION_RT_MARGIN, cfg.EXCLUSION_MZ_TOLERANCE
        )
        print(f"  exclusion list : {len(excl_entries)} entries  ->  "
              f"removed {len(removed_excl)} feature(s) for PCA")
        excl_log = os.path.join(cfg.OUTPUT_DIR, "features_removed_exclusion.csv")
        removed_excl.to_csv(excl_log, index=False)
        print(f"  -> {excl_log}")

    processed_pca, zv_pca = _transform(matrix_pca, norm_pca, log_pca, scal_pca)
    zero_var_log["pca"] = zv_pca

    out_pca = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_processed_pca.csv")
    processed_pca.to_csv(out_pca)
    print(f"  -> {out_pca}  ({processed_pca.shape[0]} samples x {processed_pca.shape[1]} features)")

    # =========================================================================
    # Zero-variance audit log
    # =========================================================================
    all_zero_var = []
    for fid in zero_var_log.get("hca", []):
        all_zero_var.append({"feature_id": fid, "analysis": "hca"})
    for fid in zero_var_log.get("pca", []):
        matched = [r for r in all_zero_var if r["feature_id"] == fid]
        if matched:
            matched[0]["analysis"] += "+pca"
        else:
            all_zero_var.append({"feature_id": fid, "analysis": "pca"})

    if all_zero_var:
        zv_df = pd.DataFrame(all_zero_var)
        meta_cols = [c for c in ("mean_rt", "mean_mz") if c in metadata.columns]
        if meta_cols:
            zv_df = (zv_df.set_index("feature_id")
                     .join(metadata[meta_cols], how="left")
                     .reset_index())
            cols  = ["feature_id"] + meta_cols + ["analysis"]
            zv_df = zv_df[[c for c in cols if c in zv_df.columns]]
        out_zv = os.path.join(cfg.OUTPUT_DIR, "features_removed_zero_variance.csv")
        zv_df.to_csv(out_zv, index=False)
        print(f"  -> {out_zv}  ({len(zv_df)} zero-variance feature(s) dropped by scaling)")

    return processed_pca, processed_hca, processed_vol


if __name__ == "__main__":
    run()
