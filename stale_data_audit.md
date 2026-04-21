# Stale-Data & CSV-Flow Audit Report
_TF_NT GC-MS Pipeline — audited 2026-04-22_

---

## Architecture Overview

`pipeline.py` creates a fresh `output/YYYY-MM-DD_HH-MM-SS/` folder, mutates
`config.OUTPUT_DIR` to that path, then calls each step's `run(config)` in
sequence.  Every step uses `os.path.join(cfg.OUTPUT_DIR, ...)` for both reads
and writes, so the timestamping is correctly propagated everywhere —
**with three concrete exceptions listed below**.

---

## Bug 1 — SILENT FEATURE-FILTER FAILURE (Execution-order race)

**File:** `blank_correction.py` lines 531–552  
**Function:** `run()`

```python
# blank_correction.py  lines 531–552
enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
if os.path.exists(enriched_path):          # ← always False in a fresh run!
    enriched = pd.read_csv(enriched_path, index_col="feature_id")
    if "molecular_formula" in enriched.columns:
        formula_map = enriched["molecular_formula"].dropna().to_dict()
```

**Why it is broken:**  
`feature_metadata_enriched.csv` is written by **Step 2c** (`compound_classification.py`).
In `pipeline.py` the execution order is:

```
Step 1  → data_import           (writes feature_metadata.csv)
Step 2  → blank_correction      (tries to read feature_metadata_enriched.csv — DOES NOT EXIST YET)
Step 2b → prevalence_histogram
Step 2c → compound_classification  (writes feature_metadata_enriched.csv — too late)
```

Because the file does not exist in the new timestamped run folder at the time
`blank_correction` runs, `os.path.exists()` returns `False`.  `formula_map` is
always `{}`.  Every keyword in `BLANK_EXCLUDE_KEYWORDS` that is meant to catch
contaminants **by molecular formula** (e.g. `"chloro"` catching `CHCl3`,
`"silan"` catching Si-containing formulas like `C10H30Si5O5`) **silently fires
on compound-name only** and misses formula matches.

This is a functional bug on every single pipeline run, with no warning printed
to the user.

### Fix — Part A: move `compound_classification` before `blank_correction` in `pipeline.py`

`compound_classification` only needs `feature_metadata.csv` (written by Step 1),
so there is no dependency conflict.

**`pipeline.py` lines 115–145 — rewrite the step order:**

```python
data_import.run(config)
print()

# Run compound classification BEFORE blank correction so that
# blank_correction can use feature_metadata_enriched.csv for
# molecular-formula keyword filtering (BLANK_EXCLUDE_KEYWORDS).
if getattr(config, "RUN_COMPOUND_CLASSIFICATION", False):
    compound_classification_step.run(config)
    print()

blank_correction.run(config)
print()
prevalence_histogram_step.run(config)
print()
if getattr(config, "RUN_CLASS_PLOTS", False):
    compound_class_plots_step.run(config)
    print()
normalization.run(config)
print()
pca_step.run(config)
print()
hca_step.run(config)
print()
if getattr(config, "RUN_HCA_DENDROGRAM", True):
    hca_dendrogram_step.run(config)
    print()
volcano_step.run(config)
print()
top_features_analysis_step.run(config)
print()
if getattr(config, "RUN_TARGETED_BOXPLOTS", True):
    targeted_boxplots_step.run(config)
    print()
blank_contaminants_report_step.run(config)
print()
classification_step.run(config)
```

### Fix — Part B: add a warning in `blank_correction.py` when the enriched file is missing

```python
# blank_correction.py  — keyword filter section  (lines 531–552)
exclude_keywords = [kw.lower() for kw in getattr(cfg, "BLANK_EXCLUDE_KEYWORDS", [])]
if exclude_keywords:
    name_map = meta["compound_name"].to_dict() \
        if meta is not None and "compound_name" in meta.columns else {}
    formula_map = {}
    enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
    if os.path.exists(enriched_path):
        enriched = pd.read_csv(enriched_path, index_col="feature_id")
        if "molecular_formula" in enriched.columns:
            formula_map = enriched["molecular_formula"].dropna().to_dict()
    else:
        print("  [warning] feature_metadata_enriched.csv not found; "
              "BLANK_EXCLUDE_KEYWORDS will match compound names only "
              "(molecular formula matching disabled).")
```

---

## Bug 2 — FILENAME MISMATCH: `compound_correlation.py` reads a file that does not exist

**File:** `compound_correlation.py` line 54  
**Function:** `run()`

```python
# compound_correlation.py  line 54
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv")   # ← wrong filename
```

**`targeted_boxplots.py` line 283** writes:

```python
out_csv = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_data.csv")   # ← actual filename
```

`compound_correlation.py` checks `os.path.exists()` and immediately returns
when the file is missing, so **it has silently produced no output on every run**
since the filename was changed.

### Fix — `compound_correlation.py` line 54

```python
# before
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv")

# after
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_data.csv")
```

---

## Bug 3 — `targeted_boxplots_unnorm.py`: top-level I/O runs at import time, bypasses timestamping

**File:** `targeted_boxplots_unnorm.py` lines 12–36 (entire script body is at module level)

```python
import config as cfg
import os

# All of this runs at MODULE LEVEL — not inside a function:
matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
peak_matrix = pd.read_csv(matrix_path)          # ← executes at import
...
targeted_df.to_csv(...)                          # ← executes at import
plt.savefig(...)                                 # ← executes at import
```

**Three problems:**

1. **Import-time side effects.** If any script ever does `import targeted_boxplots_unnorm`,
   all the file I/O and plot generation runs immediately.
2. **Bypasses timestamping.** At import time, `cfg.OUTPUT_DIR` is still `"output"`
   (the unmutated default from `config.py`).  `pipeline.py` mutates it only in
   `main()`, which runs after all module-level imports.  This script **always reads
   from and writes to the base `output/` directory**, never the timestamped run folder.
3. **Not invoked by `pipeline.py`** — its output cannot be traced to a consistent run.

### Fix — wrap all code in `run()` with a guard

```python
# targeted_boxplots_unnorm.py  — full rewrite

import math
import os
import sys

here    = os.path.dirname(os.path.abspath(__file__))
venv_py = os.path.join(here, ".venv", "bin", "python")
if os.path.exists(venv_py) and not sys.executable.startswith(os.path.join(here, ".venv")):
    os.execv(venv_py, [venv_py] + sys.argv)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

import config as _config_module


def run(cfg=_config_module):
    matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(
            f"{matrix_path} not found — run blank_correction.py first.")

    peak_matrix     = pd.read_csv(matrix_path)
    targeted_list   = cfg.EXCLUSION_LIST
    sample_cols_all = [c for c in peak_matrix.columns if c != "feature_id"]
    column_names    = ["rt", "mz", "name"] + sample_cols_all

    result_rows = []
    for target in targeted_list:
        mz_target = float(target[1])
        rt_target = float(target[0])
        row_data  = list(target[:3])
        for _, feat_row in peak_matrix.iterrows():
            fid        = feat_row["feature_id"].split("_")
            mz_feature = float(fid[0])
            rt_feature = float(fid[1])
            if (abs(rt_target - rt_feature) <= cfg.EXCLUSION_RT_MARGIN and
                    abs(mz_target  - mz_feature) <= cfg.EXCLUSION_MZ_TOLERANCE):
                row_data = list(target[:3]) + [feat_row[c] for c in sample_cols_all]
                break
        result_rows.append(row_data)

    targeted_df = pd.DataFrame(result_rows, columns=column_names)
    out_csv = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_unnorm.csv")
    targeted_df.to_csv(out_csv, index=False)

    sr_cols = [c for c in sample_cols_all if c.startswith("S-R")]
    s_cols  = [c for c in sample_cols_all
                if c.startswith("S") and not c.startswith("S-R")]

    nrows = 3
    ncols = math.ceil(len(targeted_df) / nrows)
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols * 2, nrows * 3))
    axs = np.array(axs).flatten()
    for i in range(len(targeted_df), nrows * ncols):
        axs[i].set_visible(False)

    fig.supylabel("y-axis: unnormalized area")

    for idx, row in targeted_df.iterrows():
        name = row["name"]
        SR   = np.array(row[sr_cols].values, dtype=float)
        S    = np.array(row[s_cols].values,  dtype=float)
        t_stat, p_value = ttest_ind(S, SR)
        p_signif = ("****" if p_value < 0.0001 else "***" if p_value < 0.001
                    else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns")
        ax = axs[idx]
        bp = ax.boxplot([S, SR], tick_labels=["S", "SR"], patch_artist=True)
        for patch, color in zip(bp["boxes"], ["coral", "steelblue"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for med in bp["medians"]:
            med.set_color("black")
        all_data = np.concatenate([S, SR])
        if len(all_data):
            lo, hi = np.nanmin(all_data), np.nanmax(all_data)
            pad    = (hi - lo) * 0.1 or abs(hi) * 0.1 + 0.1
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(f"{name}\n{p_signif}")
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plots_dir = os.path.join(cfg.OUTPUT_DIR, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    out_png = os.path.join(plots_dir, "targeted_boxplots_unnorm.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out_png}")


if __name__ == "__main__":
    run()
```

---

## Bug 4 — `config.py`: `PUBCHEM_CACHE_FILE` is a bare relative path

**File:** `config.py` line 532

```python
PUBCHEM_CACHE_FILE = "output/pubchem_cache.json"
```

Intentionally outside the timestamped folder (the cache persists across runs),
but as a bare relative path it breaks whenever a script is run from any directory
other than the project root.

### Fix — `config.py` line 532

```python
import os as _os
PUBCHEM_CACHE_FILE = _os.path.join(
    _os.path.dirname(_os.path.abspath(__file__)), "output", "pubchem_cache.json"
)
```

---

## Bug 5 — `classification.py`: reference files use bare relative paths

**File:** `classification.py` lines 57–60

```python
ref_path  = os.path.join("DATA", "references", "references.csv")
farn_path = os.path.join("DATA", "references", "farn1-46.csv")
ref_df  = pd.read_csv(ref_path)   # crashes if cwd != project root
farn_df = pd.read_csv(farn_path)
```

These are static reference tables, not outputs of any pipeline step, so they
cannot go stale from a previous run.  But bare relative paths will break if the
script is ever run from any other directory.

### Fix — `classification.py` lines 57–60

```python
_here     = os.path.dirname(os.path.abspath(__file__))
ref_path  = os.path.join(_here, "DATA", "references", "references.csv")
farn_path = os.path.join(_here, "DATA", "references", "farn1-46.csv")
ref_df  = pd.read_csv(ref_path)
farn_df = pd.read_csv(farn_path)
```

---

## Non-bug finding: all step-to-step data flows go through CSV round-trips

Every `run()` function returns the final DataFrame(s) it produces, but
`pipeline.py` discards every return value:

```python
# pipeline.py — returned DataFrames are thrown away
data_import.run(config)            # (no return)
blank_correction.run(config)       # returns corrected  ← discarded
normalization.run(config)          # returns (pca, hca, vol) ← discarded
pca_step.run(config)               # ...
```

Every step therefore re-reads from disk what the previous step just wrote.
This is **correct and safe** — each CSV is always written before the next step
reads it — but it is slower and more fragile (a failed write on step N causes
step N+1 to raise a `FileNotFoundError` rather than a clean in-memory error).

### Optional refactor pattern (preserves standalone-run compatibility)

Add an optional `raw=` parameter to each `run()`.  If supplied, use the
in-memory DataFrame.  If `None`, load from disk as today.

```python
# Example: blank_correction.py
def run(cfg=config, raw=None):
    ...
    if raw is not None:
        matrix = raw["peak_matrix"]
        blank_per_file = raw["blank_per_feature"]
        meta = raw.get("feature_metadata")
        group_map = raw.get("group_map", {})
    else:
        matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_raw.csv")
        matrix      = pd.read_csv(matrix_path, index_col="feature_id")
        ...
    ...
    return corrected   # already returned today

# pipeline.py
artifacts = data_import.run(config)           # returns dict of DataFrames
corrected = blank_correction.run(config, raw=artifacts)
processed = normalization.run(config, raw={"peak_matrix": corrected, ...})
...
```

---

## Summary table

| # | File | Lines | Issue | Severity |
|---|------|-------|-------|----------|
| 1 | `blank_correction.py` | 531–552 | Execution-order race: `feature_metadata_enriched.csv` never exists when blank_correction runs; molecular-formula keyword filter silently skipped on every run | **High** |
| 2 | `compound_correlation.py` | 54 | Filename mismatch: reads `targeted_boxplot.csv`, file is `targeted_boxplot_data.csv`; script silently produces no output | **High** |
| 3 | `targeted_boxplots_unnorm.py` | 12–36 | All code at module level; I/O runs at import time; reads from `output/` not the timestamped run folder | **High** |
| 4 | `config.py` | 532 | `PUBCHEM_CACHE_FILE` is a bare relative path; breaks when cwd ≠ project root | Low |
| 5 | `classification.py` | 57–60 | Reference CSV paths are bare relative paths | Low |
| 6 | `pipeline.py` | 115–145 | All DataFrames passed via CSV round-trip; opportunity for in-memory flow | Optional |

---

## Safest execution order after applying Bug 1 fix

```
1.  data_import                            → feature_metadata.csv, peak_matrix_raw.csv, ...
1b. compound_classification (if enabled)   → feature_metadata_enriched.csv   ← MOVED EARLIER
2.  blank_correction                       → peak_matrix_blank_corrected.csv  (now has enriched CSV)
2b. prevalence_histogram
2d. compound_class_plots (if enabled)      → reads enriched CSV (still available)
3.  normalization                          → 3 processed matrices
4.  pca                                    → pca_scores.csv, pca_loadings.csv
5.  hca                                    → hca heatmap
5b. hca_dendrogram (if enabled)
6.  volcano
7.  top_features_analysis
7b. targeted_boxplots (if enabled)
8.  blank_contaminants_report
9.  classification
```
