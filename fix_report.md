# Bug Fix Report
_TF_NT GC-MS Pipeline — implemented 2026-04-22_

All five bugs identified in `stale_data_audit.md` have been fixed.
Bugs are ordered high-severity first, then low-severity.

---

## Fix 1a — Execution-order race: `pipeline.py` step order

**Severity:** High  
**Files changed:** `pipeline.py`

### How it was

`compound_classification` (Step 2c) ran **after** `blank_correction` (Step 2):

```
data_import           → writes feature_metadata.csv
blank_correction      → tries to use feature_metadata_enriched.csv  ← NOT YET WRITTEN
...
compound_classification → writes feature_metadata_enriched.csv (too late)
```

### What was wrong

`blank_correction.py` contains a keyword filter (`BLANK_EXCLUDE_KEYWORDS`) that
can match contaminants by **compound name OR molecular formula**.  The molecular
formula data lives in `feature_metadata_enriched.csv`, which is produced by
`compound_classification.py`.  Because `compound_classification` ran after
`blank_correction`, the enriched file never existed in the current run folder
when `blank_correction` needed it.  `formula_map` was always empty, so
formula-based keyword matches (e.g. `"chloro"` catching `CHCl3`, `"silan"`
catching Si-containing formulas) **silently never fired on any run**.

### What was changed

`compound_classification_step.run(config)` was moved to execute immediately
after `data_import.run(config)`, before `blank_correction.run(config)`.
`compound_classification` only needs `feature_metadata.csv` (written by
`data_import`), so there is no dependency conflict.  A comment was added to
explain the ordering requirement.

**`pipeline.py` — before:**
```python
data_import.run(config)
blank_correction.run(config)
prevalence_histogram_step.run(config)
if getattr(config, "RUN_COMPOUND_CLASSIFICATION", False):
    compound_classification_step.run(config)
```

**`pipeline.py` — after:**
```python
data_import.run(config)
# compound_classification runs here (before blank_correction) so that
# feature_metadata_enriched.csv is available for molecular-formula keyword
# filtering in blank_correction (BLANK_EXCLUDE_KEYWORDS).
if getattr(config, "RUN_COMPOUND_CLASSIFICATION", False):
    compound_classification_step.run(config)
blank_correction.run(config)
prevalence_histogram_step.run(config)
```

### How it now works

When `RUN_COMPOUND_CLASSIFICATION = True`, `feature_metadata_enriched.csv` is
written into the current run folder before `blank_correction` starts.
`blank_correction` finds the file, loads the `molecular_formula` column, and
correctly filters compounds whose formula contains any keyword in
`BLANK_EXCLUDE_KEYWORDS` (in addition to the existing name matching).

When `RUN_COMPOUND_CLASSIFICATION = False`, `feature_metadata_enriched.csv` is
not created and `formula_map` remains empty — the same behaviour as before, but
now a warning is printed (see Fix 1b).

---

## Fix 1b — Execution-order race: silent failure in `blank_correction.py`

**Severity:** High  
**Files changed:** `blank_correction.py`

### How it was

```python
enriched_path = os.path.join(cfg.OUTPUT_DIR, "feature_metadata_enriched.csv")
if os.path.exists(enriched_path):
    enriched = pd.read_csv(enriched_path, ...)
    ...
# (no else branch — silent skip)
```

### What was wrong

When `feature_metadata_enriched.csv` was absent (either because
`RUN_COMPOUND_CLASSIFICATION = False` or because the old ordering placed
`compound_classification` after `blank_correction`), the code silently set
`formula_map = {}` and continued without notifying the user that half of the
keyword filter was disabled.

### What was changed

An `else` branch was added to print a clear warning whenever the enriched file
is missing so the user knows formula matching is disabled and why.

**`blank_correction.py` — after (lines 535–543):**
```python
if os.path.exists(enriched_path):
    enriched = pd.read_csv(enriched_path, index_col="feature_id")
    if "molecular_formula" in enriched.columns:
        formula_map = enriched["molecular_formula"].dropna().to_dict()
else:
    print("  [warning] feature_metadata_enriched.csv not found; "
          "BLANK_EXCLUDE_KEYWORDS will match compound names only "
          "(molecular formula matching disabled). Run compound_classification "
          "before blank_correction to enable formula matching.")
```

### How it now works

If `feature_metadata_enriched.csv` is present (the normal pipeline case with
`RUN_COMPOUND_CLASSIFICATION = True`), it is loaded silently as before.  If
absent, a one-line warning is printed and the filter continues on names only —
no silent data loss.

---

## Fix 2 — Filename mismatch: `compound_correlation.py` read a file that did not exist

**Severity:** High  
**Files changed:** `compound_correlation.py`

### How it was

```python
# compound_correlation.py  line 54 (old)
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv")
```

### What was wrong

`targeted_boxplots.py` writes its output as:

```python
out_csv = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_data.csv")
```

The filename changed at some point from `targeted_boxplot.csv` to
`targeted_boxplot_data.csv` but `compound_correlation.py` was never updated.
`os.path.exists()` returned `False` on every run, the function printed a
"not found" message and returned immediately.  **`compound_correlation.py` has
produced no output on any pipeline run since the rename.**

### What was changed

One line, correcting the filename:

**`compound_correlation.py` — before:**
```python
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot.csv")
```

**`compound_correlation.py` — after:**
```python
csv_path = os.path.join(cfg.OUTPUT_DIR, "targeted_boxplot_data.csv")
```

### How it now works

`compound_correlation.py` now finds the file written by `targeted_boxplots.py`
and correctly produces the pairwise correlation scatter grid at
`{OUTPUT_DIR}/plots/compound_correlation.png`.

---

## Fix 3 — `targeted_boxplots_unnorm.py`: module-level I/O bypassed timestamping

**Severity:** High  
**Files changed:** `targeted_boxplots_unnorm.py`

### How it was

The entire script body (file read, DataFrame build, CSV write, matplotlib figure,
`plt.savefig`) ran at **module level** — outside any function.

```python
import config as cfg
# ... (no function wrapper)
matrix_path = os.path.join(cfg.OUTPUT_DIR, "peak_matrix_blank_corrected.csv")
peak_matrix = pd.read_csv(matrix_path)          # executed at import time
...
targeted_df.to_csv(...)                          # executed at import time
plt.savefig(...)                                 # executed at import time
```

### What was wrong

Two interacting problems:

1. **Wrong `OUTPUT_DIR` at import time.**  `pipeline.py` mutates
   `config.OUTPUT_DIR` to the new timestamped run folder inside `main()`, which
   runs *after* all module-level `import` statements at the top of the file.
   So when Python imported `targeted_boxplots_unnorm` (if ever done), all I/O
   ran against `OUTPUT_DIR = "output"` — the default base directory, not the
   current run's timestamped folder.  Any output would be written to the wrong
   location, overwriting files from previous runs.

2. **Import-time side effects.**  Any future `import targeted_boxplots_unnorm`
   statement in any script would immediately read a CSV, build a plot, and write
   two files — with no way to prevent or control it.

### What was changed

The entire script was restructured into a proper `run(cfg)` function with an
`if __name__ == "__main__": run()` guard, matching the pattern used by every
other script in the pipeline.  A `FileNotFoundError` guard was added for the
input CSV.  The docstring was updated.

**`targeted_boxplots_unnorm.py` — structure before:**
```
[imports]
matrix_path = ...       # module-level I/O
peak_matrix = pd.read_csv(...)
...
plt.savefig(...)
```

**`targeted_boxplots_unnorm.py` — structure after:**
```
[imports]

def run(cfg=_config_module):
    matrix_path = os.path.join(cfg.OUTPUT_DIR, ...)   # reads correct run folder
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(...)
    peak_matrix = pd.read_csv(matrix_path)
    ...
    plt.savefig(...)
    plt.close(fig)

if __name__ == "__main__":
    run()
```

### How it now works

`run(cfg)` receives the already-mutated `config` object from `pipeline.py` (or
uses the default module when run standalone).  All file paths are resolved from
`cfg.OUTPUT_DIR` at call time, which is the correct timestamped run folder.
Importing the module no longer has any side effects.

---

## Fix 4 — `config.py`: `PUBCHEM_CACHE_FILE` was a bare relative path

**Severity:** Low  
**Files changed:** `config.py`

### How it was

```python
PUBCHEM_CACHE_FILE = "output/pubchem_cache.json"
```

### What was wrong

This is a bare relative path that is resolved relative to the **current working
directory** at runtime, not relative to the project root.  Running any pipeline
script from a different directory (e.g. `python TF_NT/pipeline.py` from the
parent folder) would cause the cache file to be created in the wrong location,
or fail to find an existing cache, forcing a full re-fetch of all PubChem data.

### What was changed

`import os as _os` was added at the top of `config.py`.  The path is now
computed at module load time using `__file__`, anchoring it to the directory
that contains `config.py` regardless of the working directory.

**`config.py` — before:**
```python
PUBCHEM_CACHE_FILE = "output/pubchem_cache.json"
```

**`config.py` — after:**
```python
import os as _os
...
PUBCHEM_CACHE_FILE = _os.path.join(
    _os.path.dirname(_os.path.abspath(__file__)), "output", "pubchem_cache.json"
)
```

### How it now works

The cache file always resolves to `<project_root>/output/pubchem_cache.json`
regardless of where the script is invoked from.  The file intentionally lives
outside timestamped run folders so it persists across pipeline runs and avoids
redundant API calls.

---

## Fix 5 — `classification.py`: reference CSVs used bare relative paths

**Severity:** Low  
**Files changed:** `classification.py`

### How it was

```python
ref_path  = os.path.join("DATA", "references", "references.csv")
farn_path = os.path.join("DATA", "references", "farn1-46.csv")
ref_df  = pd.read_csv(ref_path)
farn_df = pd.read_csv(farn_path)
```

### What was wrong

Same class of problem as Fix 4.  Bare relative paths are resolved against the
working directory.  Running `classification.py` from any directory other than
the project root would produce a `FileNotFoundError` for the reference files,
stopping the final classification step.

### What was changed

`os.path.dirname(os.path.abspath(__file__))` is stored in `_here` and used to
build absolute paths to both reference files.

**`classification.py` — before:**
```python
ref_path  = os.path.join("DATA", "references", "references.csv")
farn_path = os.path.join("DATA", "references", "farn1-46.csv")
ref_df  = pd.read_csv(ref_path)
farn_df = pd.read_csv(farn_path)
```

**`classification.py` — after:**
```python
_here     = os.path.dirname(os.path.abspath(__file__))
ref_path  = os.path.join(_here, "DATA", "references", "references.csv")
farn_path = os.path.join(_here, "DATA", "references", "farn1-46.csv")
ref_df  = pd.read_csv(ref_path)
farn_df = pd.read_csv(farn_path)
```

### How it now works

Both reference files are located relative to `classification.py` itself, so the
script works correctly regardless of the working directory from which it or
`pipeline.py` is invoked.

---

## Summary of all changes

| Bug | Severity | Files changed | Lines changed |
|-----|----------|---------------|---------------|
| 1a — Step order: `compound_classification` moved before `blank_correction` | High | `pipeline.py` | 118–129 |
| 1b — Missing warning when `feature_metadata_enriched.csv` absent | High | `blank_correction.py` | 539–543 |
| 2 — Filename mismatch in `compound_correlation.py` | High | `compound_correlation.py` | 54 |
| 3 — Module-level I/O in `targeted_boxplots_unnorm.py` | High | `targeted_boxplots_unnorm.py` | entire file |
| 4 — Bare relative path for `PUBCHEM_CACHE_FILE` | Low | `config.py` | 8–9, 532–534 |
| 5 — Bare relative paths for reference CSVs | Low | `classification.py` | 57–61 |
