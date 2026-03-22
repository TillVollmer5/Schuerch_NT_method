"""
data_import.py - Step 1 of the GCMS processing pipeline.

Loads TraceFinder CSV exports from DATA_DIR, classifies files into samples /
QC / blanks, detects features by retention-time clustering across all samples,
and writes the following files to OUTPUT_DIR:

  peak_matrix_raw.csv    - features x samples  (NaN filled with 0)
  feature_metadata.csv   - feature_id, mean_rt, mean_mz
  blank_features.csv     - feature_id, max_blank_area
  qc_matrix.csv          - features x QC samples  (saved for reference / QC checks)

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
    """
    if len(file_dict) < 2:
        return file_dict

    ref_name = sorted(file_dict.keys())[0]
    ref_df   = file_dict[ref_name]

    if "Reference m/z" not in ref_df.columns or "Retention Time" not in ref_df.columns:
        return file_dict

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

    return file_dict


# --- Feature detection --------------------------------------------------------

def _pool_peaks(file_dict, value_col):
    """
    Flatten all DataFrames in *file_dict* into a list of peak dicts.

    Each dict has keys: sample, rt, mz, area.
    Rows with unparseable numeric values are skipped.
    """
    peaks = []
    for sample_name, df in file_dict.items():
        for _, row in df.iterrows():
            try:
                peaks.append({
                    "sample": sample_name,
                    "rt":     float(row["Retention Time"]),
                    "mz":     float(row["Reference m/z"]),
                    "area":   float(row.get(value_col, 0) or 0),
                })
            except (ValueError, TypeError):
                pass
    return peaks


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
    4. Within each final cluster, if the same sample appears more than once,
       keep the peak with the highest area.

    Feature ID format
    -----------------
    use_mz=False : "RT_{mean_rt:.4f}"
    use_mz=True  : "{mean_mz:.4f}_{mean_rt:.4f}"

    Parameters
    ----------
    peaks        : list of dicts  (from _pool_peaks)
    rt_margin    : float  minutes
    use_mz       : bool
    mz_tolerance : float  Da

    Returns
    -------
    list of dicts, each with: feature_id, rt, mz, sample_areas (dict)
    """
    if not peaks:
        return []

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

    # -- build feature records -------------------------------------------------
    features = []
    for cl in rt_clusters:
        mean_rt = sum(p["rt"]  for p in cl) / len(cl)
        mean_mz = sum(p["mz"]  for p in cl) / len(cl)
        fid     = f"{mean_mz:.4f}_{mean_rt:.4f}" if use_mz else f"RT_{mean_rt:.4f}"

        # one area per sample - keep maximum when same sample appears twice
        sample_areas = {}
        for p in cl:
            s = p["sample"]
            if s not in sample_areas or p["area"] > sample_areas[s]:
                sample_areas[s] = p["area"]

        features.append({
            "feature_id":   fid,
            "rt":           mean_rt,
            "mz":           mean_mz,
            "sample_areas": sample_areas,
        })

    return features


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


def build_blank_table(features, blank_dict, rt_margin,
                      use_mz=False, mz_tolerance=0.005, value_col="Area"):
    """
    For each sample feature, find the maximum peak area across all blank files
    within +-rt_margin (and +-mz_tolerance when use_mz=True).

    Returns a Series indexed by feature_id with the maximum blank area.
    Features with no blank signal get a value of 0.
    """
    blank_peaks = _pool_peaks(blank_dict, value_col)
    max_blank   = {}

    for feat in features:
        best = 0.0
        for bp in blank_peaks:
            rt_ok = abs(bp["rt"] - feat["rt"]) <= rt_margin
            mz_ok = (not use_mz) or (abs(bp["mz"] - feat["mz"]) <= mz_tolerance)
            if rt_ok and mz_ok and bp["area"] > best:
                best = bp["area"]
        max_blank[feat["feature_id"]] = best

    s = pd.Series(max_blank, name="max_blank_area")
    s.index.name = "feature_id"
    return s


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
    if cfg.ALIGN_RT:
        samples = align_retention_times(samples, cfg.MZ_ALIGN_TOLERANCE)
        blanks  = align_retention_times(blanks,  cfg.MZ_ALIGN_TOLERANCE)
        print(f"  RT alignment    : enabled (ref = {sorted(samples.keys())[0] if samples else '-'})")

    # feature detection from all sample peaks
    sample_peaks = _pool_peaks(samples, cfg.VALUE_COL)
    features     = detect_features(
        sample_peaks,
        rt_margin    = cfg.RT_MARGIN,
        use_mz       = cfg.USE_MZ,
        mz_tolerance = cfg.MZ_TOLERANCE,
    )
    print(f"  detected {len(features)} features  "
          f"(RT margin={cfg.RT_MARGIN} min, use_mz={cfg.USE_MZ})")

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

    # feature metadata (useful for downstream annotation)
    meta = pd.DataFrame([
        {"feature_id": f["feature_id"], "mean_rt": f["rt"], "mean_mz": f["mz"]}
        for f in features
    ]).set_index("feature_id")
    out_meta = os.path.join(cfg.OUTPUT_DIR, "feature_metadata.csv")
    meta.to_csv(out_meta)
    print(f"  -> {out_meta}")

    # blank reference table
    blank_table = build_blank_table(
        features, blanks,
        rt_margin    = cfg.RT_MARGIN,
        use_mz       = cfg.USE_MZ,
        mz_tolerance = cfg.MZ_TOLERANCE,
        value_col    = cfg.VALUE_COL,
    )
    out_blank = os.path.join(cfg.OUTPUT_DIR, "blank_features.csv")
    blank_table.to_frame().to_csv(out_blank)
    print(f"  -> {out_blank}  "
          f"({(blank_table > 0).sum()} features with blank signal)")

    return matrix, blank_table, group_map


if __name__ == "__main__":
    run()
