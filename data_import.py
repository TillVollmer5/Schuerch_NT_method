"""
data_import.py - Step 1 of the GCMS processing pipeline.

Loads TraceFinder CSV exports from DATA_DIR, classifies files into samples /
QC / blanks, detects features by retention-time clustering across all samples,
and writes the following files to OUTPUT_DIR:

  peak_matrix_raw.csv        - features x samples  (missing values filled with 0)
  feature_metadata.csv       - feature_id, mean_rt, mean_mz, cluster spread stats,
                               n_samples_detected, n_contributing_peaks
  blank_features.csv         - feature_id, max_blank_area
  feature_peak_log.csv       - full provenance: one row per raw peak,
                               which feature it was assigned to, source sample,
                               ref_mz, rt_raw (pre-alignment), rt_aligned,
                               rt_shift applied, area, and whether the peak was
                               selected (True) or replaced by a higher-area
                               duplicate from the same sample in the same cluster
  rt_alignment_shifts.csv    - median RT shift applied per sample (0.0 if none)

Usage:
    python data_import.py
"""

import os
import sys

# re-exec under workspace virtualenv if not already active
here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import glob
import pandas as pd

import config


# --- File loading -------------------------------------------------------------

def load_files(data_dir, blank_prefix, sample_groups):
    """
    Load and classify all CSV files in *data_dir*.

    Classification is driven by config parameters:
      - Files whose name starts with *blank_prefix* -> blanks (excluded from matrix)
      - All other files are matched against *sample_groups* (ordered list of
        (group_name, prefix) pairs); the first matching prefix wins.
      - Unmatched files are placed in group "unknown".

    Returns
    -------
    samples   : dict  stem -> DataFrame
    blanks    : dict  stem -> DataFrame
    group_map : dict  stem -> group_name
    """
    samples, blanks, group_map = {}, {}, {}

    for fp in sorted(glob.glob(os.path.join(data_dir, "*.csv"))):
        name = os.path.basename(fp)
        stem = os.path.splitext(name)[0]

        df = _read_csv(fp)
        if df is None:
            print(f"  [skip] could not parse {name}")
            continue

        if name.lower().startswith(blank_prefix.lower()):
            blanks[stem] = df
            continue

        # match against sample groups in order (most specific prefix first)
        group = "unknown"
        for group_name, prefix in sample_groups:
            if name.upper().startswith(prefix.upper()):
                group = group_name
                break

        samples[stem] = df
        group_map[stem] = group

    return samples, blanks, group_map


def _read_csv(fp):
    """
    Read a TraceFinder CSV, handling an optional non-data header row.
    Returns a DataFrame on success, None on failure.
    """
    required = {"Reference m/z", "Retention Time"}
    for skiprows in (1, 0):
        try:
            df = pd.read_csv(fp, skiprows=skiprows, engine="python", on_bad_lines="warn")
            if required.issubset(df.columns):
                return df
        except Exception:
            pass
    return None


# --- RT alignment (optional pre-processing) -----------------------------------

def align_retention_times(file_dict, mz_tolerance=0.1):
    """
    Apply a median RT-shift correction so that the same compound in different
    samples lands at the same retention time.

    The first file (alphabetically) is used as the alignment reference.
    Alignment is skipped silently for files missing the required columns.

    Returns
    -------
    file_dict : dict  (modified in-place with corrected RTs)
    shifts    : dict  sample_name -> float  (median RT shift applied; 0.0 if none)
    """
    shifts = {name: 0.0 for name in file_dict}

    if len(file_dict) < 2:
        return file_dict, shifts

    ref_name = sorted(file_dict.keys())[0]
    ref_df   = file_dict[ref_name]

    if "Reference m/z" not in ref_df.columns or "Retention Time" not in ref_df.columns:
        return file_dict, shifts

    ref_sorted = ref_df.sort_values("Reference m/z").reset_index(drop=True)

    for name, df in file_dict.items():
        if name == ref_name:
            continue
        if "Reference m/z" not in df.columns or "Retention Time" not in df.columns:
            continue

        sample_sorted = df.sort_values("Reference m/z").reset_index(drop=True)
        merged = pd.merge_asof(
            sample_sorted,
            ref_sorted,
            on="Reference m/z",
            tolerance=mz_tolerance,
            direction="nearest",
            suffixes=("", "_ref"),
        )
        if "Retention Time_ref" in merged.columns:
            shift = (merged["Retention Time_ref"] - merged["Retention Time"]).median()
            if pd.notna(shift):
                file_dict[name] = df.copy()
                file_dict[name]["Retention Time"] = df["Retention Time"] + shift
                shifts[name] = float(shift)

    return file_dict, shifts


# --- Feature detection --------------------------------------------------------

def _pool_peaks(file_dict, value_col, rt_shifts=None, name_col=None):
    """
    Flatten all DataFrames in *file_dict* into a list of peak dicts.

    Each dict has keys: sample, rt (aligned), rt_raw (pre-alignment),
    rt_aligned, rt_shift, mz, area, name.
    Rows with unparseable numeric values are skipped.

    Parameters
    ----------
    file_dict  : dict  sample_name -> DataFrame
    value_col  : str   column to extract for area
    rt_shifts  : dict  sample_name -> float  (from align_retention_times)
                       If None, all shifts are treated as 0.
    name_col   : str or None  column to extract compound name from
                       (e.g. "Name"); None or "" = skip name extraction
    """
    peaks = []
    for sample_name, df in file_dict.items():
        shift = (rt_shifts or {}).get(sample_name, 0.0)
        for _, row in df.iterrows():
            try:
                rt_aligned = float(row["Retention Time"])
                name_val = ""
                if name_col:
                    try:
                        name_val = str(row.get(name_col, "") or "").strip()
                    except Exception:
                        pass
                peaks.append({
                    "sample":     sample_name,
                    "rt":         rt_aligned,          # used for clustering
                    "rt_raw":     rt_aligned - shift,  # original RT before alignment
                    "rt_aligned": rt_aligned,
                    "rt_shift":   shift,
                    "mz":         float(row["Reference m/z"]),
                    "area":       float(row.get(value_col, 0) or 0),
                    "name":       name_val,
                })
            except (ValueError, TypeError):
                pass
    return peaks


def _split_on_conflict(cluster):
    """
    Recursively split *cluster* until no sub-cluster contains two peaks from
    the same sample.

    Rationale: a single sample injection can only detect each compound once.
    If two peaks from the same sample fall in the same RT cluster, they must
    represent two chemically distinct features — the cluster must be split.

    Split criterion: the largest RT gap between consecutive peaks (sorted by
    RT).  Splitting at the widest gap keeps peaks that are closest in RT
    together, which is consistent with the RT-similarity assumption that
    underlies the original greedy clustering.

    Recursion terminates when every sub-cluster has unique sample membership,
    or when a sub-cluster contains only one peak (cannot be split further).

    Returns a list of sub-clusters (each is a list of peak dicts).
    """
    samples = [p["sample"] for p in cluster]
    if len(samples) == len(set(samples)):
        return [cluster]          # no conflict — nothing to do

    if len(cluster) <= 1:
        return [cluster]          # degenerate — cannot split further

    sorted_cl = sorted(cluster, key=lambda p: p["rt"])
    gaps      = [sorted_cl[i + 1]["rt"] - sorted_cl[i]["rt"]
                 for i in range(len(sorted_cl) - 1)]
    split_at  = gaps.index(max(gaps))   # index of peak just before the gap

    left  = sorted_cl[:split_at + 1]
    right = sorted_cl[split_at + 1:]

    return _split_on_conflict(left) + _split_on_conflict(right)


def detect_features(peaks, rt_margin, use_mz=False, mz_tolerance=0.005):
    """
    Cluster peaks from multiple samples into features.

    Algorithm
    ---------
    1. Sort all peaks by RT.
    2. Greedy RT clustering: start a new cluster when the current peak's RT
       exceeds the cluster's anchor RT by more than *rt_margin*.
    3. Optional m/z sub-clustering: split each RT cluster by m/z proximity
       (within *mz_tolerance*) to separate co-eluting compounds.
    4. Same-sample conflict resolution: if a cluster contains two or more
       peaks from the same sample, it must contain at least two distinct
       features.  Split recursively at the largest RT gap until every
       sub-cluster has unique sample membership.
    5. Build one feature record per final cluster.  Every peak in each
       cluster is marked selected=True (no data is discarded).

    Feature ID format
    -----------------
    use_mz=False : "RT_{mean_rt:.4f}"
    use_mz=True  : "{mean_mz:.4f}_{mean_rt:.4f}"

    Each feature record includes a *peak_log* list with one entry per raw
    peak assigned to this cluster.  After conflict resolution every peak is
    selected=True (each sample appears at most once per cluster).

    Returns
    -------
    features    : list of dicts, each with:
                      feature_id, rt, mz, compound_name, sample_areas (dict),
                      peak_log (list of dicts), rt_values (list), mz_values (list)
    n_splits    : int  number of clusters that were split by conflict resolution
    """
    if not peaks:
        return [], 0

    sorted_peaks = sorted(peaks, key=lambda x: x["rt"])

    # -- greedy RT clustering --------------------------------------------------
    rt_clusters, cluster = [], [sorted_peaks[0]]
    for peak in sorted_peaks[1:]:
        if peak["rt"] - cluster[0]["rt"] <= rt_margin:
            cluster.append(peak)
        else:
            rt_clusters.append(cluster)
            cluster = [peak]
    rt_clusters.append(cluster)

    # -- optional m/z sub-clustering ------------------------------------------
    if use_mz:
        sub_clusters = []
        for cl in rt_clusters:
            mz_sorted = sorted(cl, key=lambda x: x["mz"])
            sub = [mz_sorted[0]]
            for p in mz_sorted[1:]:
                if p["mz"] - sub[0]["mz"] <= mz_tolerance:
                    sub.append(p)
                else:
                    sub_clusters.append(sub)
                    sub = [p]
            sub_clusters.append(sub)
        rt_clusters = sub_clusters

    # -- same-sample conflict resolution --------------------------------------
    # Two peaks from the same sample cannot be the same feature.  Any cluster
    # containing a duplicate sample is split at the largest RT gap, recursively,
    # until all sub-clusters have unique sample membership.
    n_splits    = 0
    resolved    = []
    for cl in rt_clusters:
        sub = _split_on_conflict(cl)
        if len(sub) > 1:
            n_splits += len(sub) - 1
        resolved.extend(sub)
    rt_clusters = resolved

    # -- build feature records -------------------------------------------------
    features = []
    for cl in rt_clusters:
        mean_rt = sum(p["rt"]  for p in cl) / len(cl)
        mean_mz = sum(p["mz"]  for p in cl) / len(cl)
        fid     = f"{mean_mz:.4f}_{mean_rt:.4f}" if use_mz else f"RT_{mean_rt:.4f}"

        # one peak per sample (guaranteed by conflict resolution)
        sample_areas = {p["sample"]: p["area"] for p in cl}

        # compound name from the peak with the highest area across all samples
        best_peak     = max(cl, key=lambda p: p.get("area", 0))
        compound_name = best_peak.get("name", "")

        # peak log — every peak is selected (no within-cluster duplicates remain)
        peak_log = [{
            "feature_id": fid,
            "sample":     p["sample"],
            "ref_mz":     p["mz"],
            "rt_raw":     p.get("rt_raw",     p["rt"]),
            "rt_aligned": p.get("rt_aligned", p["rt"]),
            "rt_shift":   p.get("rt_shift",   0.0),
            "area":       p["area"],
            "selected":   True,
        } for p in cl]

        features.append({
            "feature_id":    fid,
            "rt":            mean_rt,
            "mz":            mean_mz,
            "compound_name": compound_name,
            "sample_areas":  sample_areas,
            "peak_log":      peak_log,
            "rt_values":     [p["rt"]  for p in cl],
            "mz_values":     [p["mz"]  for p in cl],
        })

    return features, n_splits


# --- Matrix and table construction -------------------------------------------

def build_matrix(features, sample_names):
    """
    Build a features x samples DataFrame from detected features.
    Missing values (feature not detected in a sample) are filled with 0.
    """
    rows   = {f["feature_id"]: f["sample_areas"] for f in features}
    matrix = pd.DataFrame(rows).T
    matrix = matrix.reindex(columns=sample_names)
    matrix.fillna(0, inplace=True)
    matrix.index.name = "feature_id"
    return matrix


def build_blank_table(features, blank_dict, rt_margin, value_col="Area"):
    """
    For each sample feature, find the blank peak with the highest area
    within +-rt_margin (RT-only matching).  The matched blank peak's
    m/z and RT are also stored so that blank_correction.py can optionally
    apply an additional m/z proximity gate (BLANK_USE_MZ in config.py).

    RT-only matching is intentional here: the m/z gate is a separate,
    independently configurable step in blank_correction.py.  This avoids
    conflating feature-detection m/z clustering (USE_MZ) with blank-matching
    m/z stringency (BLANK_USE_MZ).

    Returns
    -------
    DataFrame indexed by feature_id with columns:
        max_blank_area  - highest blank area within the RT window (0 if no match)
        blank_rt        - RT of the matched blank peak (NaN if no match)
        blank_mz        - m/z of the matched blank peak (NaN if no match)
    """
    blank_peaks = _pool_peaks(blank_dict, value_col)   # shifts default to 0
    rows = {}

    for feat in features:
        best_area = 0.0
        best_mz   = float("nan")
        best_rt   = float("nan")
        for bp in blank_peaks:
            if abs(bp["rt"] - feat["rt"]) <= rt_margin and bp["area"] > best_area:
                best_area = bp["area"]
                best_mz   = bp["mz"]
                best_rt   = bp["rt"]
        rows[feat["feature_id"]] = {
            "max_blank_area": best_area,
            "blank_rt":       best_rt,
            "blank_mz":       best_mz,
        }

    df = pd.DataFrame(rows).T
    df.index.name = "feature_id"
    return df


# --- Main ---------------------------------------------------------------------

def run(cfg=config):
    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)

    print("-- Step 1: data import -------------------------------------------")

    # load and classify files using config-defined groups
    samples, blanks, group_map = load_files(
        cfg.DATA_DIR, cfg.BLANK_PREFIX, cfg.SAMPLE_GROUPS
    )

    # print group summary
    for group_name, _ in cfg.SAMPLE_GROUPS:
        members = sorted(s for s, g in group_map.items() if g == group_name)
        print(f"  group '{group_name}': {members}")
    unknown = sorted(s for s, g in group_map.items() if g == "unknown")
    if unknown:
        print(f"  group 'unknown': {unknown}")
    print(f"  blanks          : {sorted(blanks.keys())}")

    # optional RT alignment (reference = first sample alphabetically)
    sample_shifts = {s: 0.0 for s in samples}
    blank_shifts  = {b: 0.0 for b in blanks}

    if cfg.ALIGN_RT:
        samples, sample_shifts = align_retention_times(samples, cfg.MZ_ALIGN_TOLERANCE)
        blanks,  blank_shifts  = align_retention_times(blanks,  cfg.MZ_ALIGN_TOLERANCE)
        ref = sorted(samples.keys())[0] if samples else "-"
        print(f"  RT alignment    : enabled (ref = {ref})")

        # save alignment shifts for traceability
        all_shifts = {**sample_shifts, **blank_shifts}
        shifts_df  = pd.DataFrame(
            [{"sample": s, "rt_shift_applied": v} for s, v in sorted(all_shifts.items())]
        )
        out_shifts = os.path.join(cfg.OUTPUT_DIR, "rt_alignment_shifts.csv")
        shifts_df.to_csv(out_shifts, index=False)
        print(f"  -> {out_shifts}")
    else:
        # write a zero-shifts file so downstream tools always find the file
        all_shifts = {**sample_shifts, **blank_shifts}
        shifts_df  = pd.DataFrame(
            [{"sample": s, "rt_shift_applied": 0.0} for s in sorted(all_shifts)]
        )
        out_shifts = os.path.join(cfg.OUTPUT_DIR, "rt_alignment_shifts.csv")
        shifts_df.to_csv(out_shifts, index=False)

    # feature detection from all sample peaks (pass shifts for RT provenance)
    name_col     = getattr(cfg, "COMPOUND_NAME_COL", "Name") or ""
    sample_peaks = _pool_peaks(samples, cfg.VALUE_COL, rt_shifts=sample_shifts,
                               name_col=name_col)
    features, n_splits = detect_features(
        sample_peaks,
        rt_margin    = cfg.RT_MARGIN,
        use_mz       = cfg.USE_MZ,
        mz_tolerance = cfg.MZ_TOLERANCE,
    )
    print(f"  detected {len(features)} features  "
          f"(RT margin={cfg.RT_MARGIN} min, use_mz={cfg.USE_MZ})")
    if n_splits > 0:
        print(f"  same-sample conflict splits : {n_splits} cluster(s) split "
              f"(peaks from the same sample separated into distinct features)")

    # peak matrix (features x all samples, ordered by group then name)
    sample_order = sorted(
        samples.keys(),
        key=lambda s: (
            [g for g, _ in cfg.SAMPLE_GROUPS].index(group_map.get(s, "unknown"))
            if group_map.get(s, "unknown") in [g for g, _ in cfg.SAMPLE_GROUPS]
            else len(cfg.SAMPLE_GROUPS),
            s,
        )
    )
    matrix     = build_matrix(features, sample_order)
    out_matrix = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
    matrix.to_csv(out_matrix)
    print(f"  -> {out_matrix}  "
          f"({matrix.shape[0]} features x {matrix.shape[1]} samples)")

    # group labels - saved alongside the matrix for use in PCA colouring etc.
    group_df = pd.DataFrame(
        {"sample": sample_order, "group": [group_map[s] for s in sample_order]}
    ).set_index("sample")
    out_groups = os.path.join(cfg.OUTPUT_DIR, "sample_groups.csv")
    group_df.to_csv(out_groups)
    print(f"  -> {out_groups}")

    # --- feature metadata (enhanced with cluster spread and detection counts) --
    n_samples = len(sample_order)
    meta_rows = []
    for f in features:
        rt_vals = f["rt_values"]
        mz_vals = f["mz_values"]
        n_det   = int((matrix.loc[f["feature_id"]] > 0).sum())
        meta_rows.append({
            "feature_id":          f["feature_id"],
            "compound_name":       f.get("compound_name", ""),
            "mean_rt":             f["rt"],
            "mean_mz":             f["mz"],
            "rt_min":              min(rt_vals),
            "rt_max":              max(rt_vals),
            "rt_std":              (sum((v - f["rt"]) ** 2 for v in rt_vals) / len(rt_vals)) ** 0.5
                                   if len(rt_vals) > 1 else 0.0,
            "mz_min":              min(mz_vals),
            "mz_max":              max(mz_vals),
            "mz_std":              (sum((v - f["mz"]) ** 2 for v in mz_vals) / len(mz_vals)) ** 0.5
                                   if len(mz_vals) > 1 else 0.0,
            "n_samples_detected":  n_det,
            "n_samples_total":     n_samples,
            "n_contributing_peaks": len(f["peak_log"]),
        })
    meta = pd.DataFrame(meta_rows).set_index("feature_id")
    out_meta = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")
    meta.to_csv(out_meta)
    print(f"  -> {out_meta}  (incl. cluster spread and detection counts)")

    # compound name map - feature_id -> compound_name (for plot labelling)
    name_map = meta[["compound_name"]].copy()
    out_name_map = os.path.join(cfg.OUTPUT_DIR, "feature_name_map.csv")
    name_map.to_csv(out_name_map)
    n_named = int((name_map["compound_name"].str.strip() != "").sum())
    print(f"  -> {out_name_map}  ({n_named} features with compound names)")

    # --- peak provenance log --------------------------------------------------
    peak_log_rows = []
    for f in features:
        peak_log_rows.extend(f["peak_log"])

    peak_log_df = pd.DataFrame(peak_log_rows, columns=[
        "feature_id", "sample", "ref_mz",
        "rt_raw", "rt_aligned", "rt_shift", "area", "selected",
    ])
    out_peak_log = os.path.join(cfg.OUTPUT_DIR, "feature_peak_log.csv")
    peak_log_df.to_csv(out_peak_log, index=False)
    n_duplicates = int((~peak_log_df["selected"]).sum())
    print(f"  -> {out_peak_log}  "
          f"({len(peak_log_df)} peaks total, {n_duplicates} within-cluster duplicates)")

    # blank reference table (RT-only matching; m/z gate is in blank_correction.py)
    blank_table = build_blank_table(
        features, blanks,
        rt_margin = cfg.RT_MARGIN,
        value_col = cfg.VALUE_COL,
    )
    out_blank = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")
    blank_table.to_csv(out_blank)
    print(f"  -> {out_blank}  "
          f"({(blank_table['max_blank_area'] > 0).sum()} features with blank signal)")

    return matrix, blank_table, group_map


if __name__ == "__main__":
    run()
