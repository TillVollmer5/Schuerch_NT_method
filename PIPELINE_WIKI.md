# Pipeline Wiki — Nontargeted GC-MS Metabolomics (TF_NT)

A detailed, step-by-step description of every statistical, mathematical, and
mass-spectrometry-relevant process in the pipeline. Written for the scientific
reader who wants to understand *what* happens to data between TraceFinder export
and the final plots.

---

## Table of Contents

1. [Overview and data flow](#1-overview-and-data-flow)
2. [Step 1 — Data import and feature detection](#2-step-1--data-import-and-feature-detection)
   - 2.1 [File loading and sample classification](#21-file-loading-and-sample-classification)
   - 2.2 [Retention-time alignment](#22-retention-time-alignment)
   - 2.3 [Peak pooling](#23-peak-pooling)
   - 2.4 [Feature detection by RT clustering](#24-feature-detection-by-rt-clustering)
   - 2.5 [Optional m/z sub-clustering](#25-optional-mz-sub-clustering)
   - 2.6 [Within-cluster duplicate handling](#26-within-cluster-duplicate-handling)
   - 2.7 [Feature ID and metadata](#27-feature-id-and-metadata)
   - 2.8 [Blank reference table](#28-blank-reference-table)
   - 2.9 [Peak matrix construction](#29-peak-matrix-construction)
3. [Step 2 — Blank correction](#3-step-2--blank-correction)
   - 3.1 [Optional m/z gate](#31-optional-mz-gate)
   - 3.2 [Fold-change filter](#32-fold-change-filter)
4. [Step 2b — Prevalence histogram](#4-step-2b--prevalence-histogram)
   - 4.1 [Prevalence calculation](#41-prevalence-calculation)
   - 4.2 [Bin alignment and x-axis scaling](#42-bin-alignment-and-x-axis-scaling)
   - 4.3 [Threshold reference lines](#43-threshold-reference-lines)
   - 4.4 [Outputs](#44-outputs)
5. [Step 3 — Normalization, transformation, and scaling](#5-step-3--normalization-transformation-and-scaling)
   - 4.1 [Prevalence filter](#41-prevalence-filter)
   - 4.2 [Exclusion list (PCA matrix only)](#42-exclusion-list-pca-matrix-only)
   - 4.3 [Sum normalization](#43-sum-normalization)
   - 4.4 [Median normalization (alternative)](#44-median-normalization-alternative)
   - 4.5 [Log transformation](#45-log-transformation)
   - 4.6 [Pareto scaling (recommended)](#46-pareto-scaling-recommended)
   - 4.7 [Auto scaling (alternative)](#47-auto-scaling-alternative)
   - 4.8 [Two output matrices](#48-two-output-matrices)
6. [Step 4 — Principal Component Analysis (PCA)](#6-step-4--principal-component-analysis-pca)
   - 6.1 [Algorithm](#61-algorithm)
   - 6.2 [Scores plot and confidence ellipses](#62-scores-plot-and-confidence-ellipses)
   - 6.3 [Loadings plot and bar chart](#63-loadings-plot-and-bar-chart)
7. [Step 5 — Hierarchical Cluster Analysis (HCA)](#7-step-5--hierarchical-cluster-analysis-hca)
   - 7.1 [Algorithm and linkage](#71-algorithm-and-linkage)
   - 7.2 [Heatmap construction](#72-heatmap-construction)
8. [Step 6 — Volcano plot](#8-step-6--volcano-plot)
   - 8.1 [Normalization for volcano (no Pareto scaling)](#81-normalization-for-volcano-no-pareto-scaling)
   - 8.2 [Log2 fold change](#82-log2-fold-change)
   - 8.3 [Mann-Whitney U test](#83-mann-whitney-u-test)
   - 8.4 [Benjamini-Hochberg FDR correction](#84-benjamini-hochberg-fdr-correction)
   - 8.5 [Classification and thresholds](#85-classification-and-thresholds)
9. [Report — Blank contaminants](#9-report--blank-contaminants)
10. [Report — Top PCA features](#10-report--top-pca-features)
11. [Feature labelling (compound names)](#11-feature-labelling-compound-names)
12. [Assumptions summary](#12-assumptions-summary)
13. [Key parameter reference](#13-key-parameter-reference)

---

## 1. Overview and data flow

```
TraceFinder CSV exports (one file per injection)
        |
        v
Step 1  data_import.py
        - classify files (sample / blank)
        - optional RT alignment
        - feature detection by RT clustering
        - build peak matrix (features x samples)
        - build blank reference table
        |
        v  peak_matrix_raw.csv
           blank_features.csv         (max blank area per feature, backward-compat)
           blank_per_feature.csv      (per-blank-file signal; enables mean/each modes)
           feature_metadata.csv
           feature_name_map.csv
           feature_peak_log.csv
           rt_alignment_shifts.csv
           sample_groups.csv
        |
        v
Step 2  blank_correction.py
        - optional m/z gate on blank matching (per blank file)
        - fold-change filter (configurable sample and blank reference modes)
        - produces full per-comparison audit log
        |
        v  peak_matrix_blank_corrected.csv
           blank_correction_audit.csv
           features_removed_blank.csv (backward-compat summary)
        |
        v
Step 2b prevalence_histogram.py
        - reads blank-corrected matrix
        - computes per-feature detection fraction across all samples
        - produces histogram with bin edges aligned to 1/N_samples steps
        - annotates PCA / HCA prevalence threshold lines from config.py
        -> prevalence_histogram.png
        -> prevalence_summary.csv
        |
        v
Step 3  normalization.py
        - prevalence filter (separately for HCA and PCA)
        - exclusion list (PCA matrix only)
        - normalization (sum or median)
        - log2(x + 1) transformation
        - Pareto or auto scaling
        - produces TWO output matrices
        |
        v  peak_matrix_processed.csv      (full feature set -> HCA, volcano)
           peak_matrix_processed_pca.csv  (filtered         -> PCA)
        |
        +--------> Step 4  pca.py        -> scores/loadings plots
        |
        +--------> Step 5  hca.py        -> clustered heatmap
        |
        +--------> Step 6  volcano.py    -> per-comparison volcano plots
                   (reads blank_corrected matrix directly, re-normalises internally)
        |
        v
Step 7  blank_contaminants_report.py
        - reads blank_correction_audit.csv, feature_name_map.csv
        - per-(feature, sample_unit, blank) traceability report:
          which blank triggered each removal, fold change, m/z delta, etc.
        -> blank_contaminants_report.csv
        |
        v
Step 8  top_features_analysis.py
        - reads pca_loadings.csv, peak_matrix_raw.csv, feature_metadata.csv,
          and raw TraceFinder CSVs
        - extracts top N features by PCA loading magnitude with compound
          names looked up from raw data by RT matching
        -> top_features_analysis.csv
```

The pipeline is deliberately linear through Steps 1-3 (each step depends on
the previous), then parallel at Steps 4-6 (independent analyses sharing the
same processed matrices). Steps 7 and 8 are post-analysis reporting steps
that depend on the outputs of earlier steps but do not modify any shared files.

---

## 2. Step 1 — Data import and feature detection

**Script:** `data_import.py`
**Input:** CSV files in `DATA/`
**Output:** `peak_matrix_raw.csv`, `feature_metadata.csv`, `feature_name_map.csv`,
`feature_peak_log.csv`, `rt_alignment_shifts.csv`, `blank_features.csv`,
`sample_groups.csv`

### 2.1 File loading and sample classification

Every `.csv` file in `DATA_DIR` is read and classified in three stages:

1. **Header detection.** TraceFinder exports sometimes prepend a non-data
   metadata row before the column headers. The reader tries `skiprows=1` first,
   then `skiprows=0`, accepting whichever parse yields a DataFrame containing
   both `"Reference m/z"` and `"Retention Time"` columns. Files that cannot be
   parsed this way are skipped.

2. **Blank classification.** Files whose name begins with `BLANK_PREFIX`
   (case-insensitive match) are assigned to the blank pool. Blanks are not
   included in the sample matrix but are used exclusively for blank correction
   in Step 2.

3. **Sample group assignment.** All remaining files are matched against
   `SAMPLE_GROUPS` — an ordered list of `(group_name, prefix)` pairs. The
   first matching prefix (case-insensitive) determines the group. Files that
   match no prefix are placed in group `"unknown"`. The ordering matters: more
   specific prefixes (e.g. `"S-R"`) must appear before broader ones (e.g.
   `"S"`) to prevent the broader prefix from consuming both groups.

**Assumption:** The file naming convention uniquely identifies biological groups.
Files named identically except for a trailing digit (e.g. `S-R1.csv`, `S-R2.csv`,
`S-R3.csv`) are replicates of the same group.

---

### 2.2 Retention-time alignment

**Config:** `ALIGN_RT`, `MZ_ALIGN_TOLERANCE`

GC-MS instruments experience session-to-session and run-to-run retention-time
drift. Even within a single sequence, early injections often elute slightly
differently from late ones due to column bleed, temperature fluctuations, or
carrier-gas pressure changes. Without correction, the same compound appearing
at (for example) 5.001 min in sample 1 and 5.043 min in sample 6 may be
assigned to two separate features.

**Algorithm (when `ALIGN_RT = True`):**

1. The file whose name sorts first alphabetically is chosen as the alignment
   reference (reference sample).
2. For every other sample, its peaks are matched to the reference peaks by m/z
   using `pandas.merge_asof` (nearest-neighbour merge on `"Reference m/z"` with
   `MZ_ALIGN_TOLERANCE` Da tolerance). This gives pairs of
   (sample RT, reference RT) for peaks that share the same nominal m/z.
3. The **median** of `(reference RT − sample RT)` across all matched peak pairs
   is computed. This median shift is robust to outliers caused by mismatched
   peaks or compound-specific RT differences.
4. The shift is added to every peak RT in that sample:
   `RT_aligned = RT_raw + shift`
5. All downstream feature detection uses the aligned RTs. The raw RT and the
   applied shift are both stored in the peak log for traceability.

**When `ALIGN_RT = False`:** all shifts are set to `0.0` and the raw RTs are
used as-is. A zero-shifts file is still written so that downstream audit tools
always find `rt_alignment_shifts.csv`.

**Assumptions:**
- RT drift is approximately constant across the entire chromatogram within a
  single injection (a single shift value per sample suffices).
- There are enough shared peaks between the reference and each sample within
  the m/z tolerance to compute a reliable median.
- The reference sample itself has zero drift (it is the reference by definition).

---

### 2.3 Peak pooling

After alignment, all sample DataFrames are flattened into a single list of peak
dictionaries. Each entry carries:

| Field | Description |
|-------|-------------|
| `sample` | Source file stem (e.g. `"S-R1"`) |
| `rt` | Aligned RT (minutes) — used for all clustering |
| `rt_raw` | Pre-alignment RT = `rt_aligned − shift` |
| `rt_aligned` | Same as `rt` (explicit alias for the provenance log) |
| `rt_shift` | Median shift applied to this sample |
| `mz` | Reference m/z (Da) |
| `area` | Chromatographic peak area (or height if `VALUE_COL = "Height"`) |
| `name` | Compound name from the `COMPOUND_NAME_COL` column (`""` if absent) |

Rows with non-numeric RT, m/z, or area values are silently skipped.

**On Area vs Height:**
Area (the default) integrates signal over time and is less sensitive to peak
shape distortions caused by co-elution or tailing. Height is faster to compute
but depends on chromatographic peak shape being consistent across runs. Area is
the standard choice for quantitative untargeted GC-MS.

---

### 2.4 Feature detection by RT clustering

Feature detection addresses the fundamental problem of defining what constitutes
"the same metabolite signal" across multiple samples. Because TraceFinder
reports one row per detected peak per file, two peaks in two samples that
represent the same compound must be identified and grouped together.

**Algorithm — greedy RT clustering:**

1. All peaks from all samples are sorted by ascending RT.
2. The first peak starts cluster 1. Its RT is the cluster **anchor**.
3. For each subsequent peak (in RT order):
   - If `peak.rt − cluster_anchor.rt ≤ RT_MARGIN` → append to the current cluster.
   - Otherwise → close the current cluster, start a new one with this peak as anchor.
4. The process produces a set of RT clusters, each spanning at most `RT_MARGIN`
   minutes from its first peak.

The greedy approach is O(n) after sorting and handles variable cluster sizes
naturally. The anchor is always the **first peak** in sorted RT order, not the
mean, which means the cluster boundary is `[anchor_rt, anchor_rt + RT_MARGIN]`.

**Key design choices and limitations:**
- `RT_MARGIN` is the only free parameter. Too small → the same compound in
  samples with slightly different RTs ends up in separate features (false
  splitting). Too large → two nearby co-eluting compounds are merged into one
  feature (false fusion).
- The greedy anchor approach does not re-centre the cluster. If peaks
  accumulate at the boundary, a cluster's effective width can be slightly less
  than `2 × RT_MARGIN` but never more than `RT_MARGIN`.
- Peaks are pooled *across all samples*, so a cluster can contain multiple
  peaks from the same sample (see §2.6).

**Assumptions:**
- The same compound elutes within `RT_MARGIN` minutes across all samples,
  after RT alignment.
- Compounds eluting closer than `RT_MARGIN` to each other are considered the
  same feature unless `USE_MZ = True` (§2.5).
- A feature may be absent from some samples (area = 0 in the matrix), which
  represents a true biological absence or a signal below the detection limit.

---

### 2.5 Optional m/z sub-clustering

**Config:** `USE_MZ`, `MZ_TOLERANCE`

In complex matrices (e.g. biological extracts with dense chromatograms), two
chemically distinct compounds may elute within `RT_MARGIN` of each other. When
`USE_MZ = True`, each RT cluster is further divided by m/z proximity:

1. Peaks within the RT cluster are sorted by ascending m/z.
2. A second greedy clustering step runs on m/z: a new sub-cluster starts when
   `peak.mz − sub_cluster_anchor.mz > MZ_TOLERANCE`.
3. Each sub-cluster becomes an independent feature.

**On Reference m/z in GC-MS:**
In TraceFinder GC-MS workflows, `"Reference m/z"` is the nominal mass (typically
the base peak or a quantitation ion) selected during compound library matching
or deconvolution. For Orbitrap data, this value has sufficient mass accuracy
(sub-5 ppm) to distinguish compounds with different elemental compositions.
`MZ_TOLERANCE = 0.005 Da` is appropriate for Orbitrap-acquired data.

**Assumption (when `USE_MZ = False`):** A single RT window is sufficient to
separate all compounds of interest, or the user is not concerned about
co-eluting compounds being merged.

---

### 2.6 Same-sample conflict resolution

If two peaks from the **same sample** fall within the same RT cluster, they
cannot represent the same feature — a single injection produces each compound
exactly once. The cluster must therefore contain at least two chemically distinct
features and must be split.

**Algorithm — recursive largest-gap split:**

1. Check whether any sample appears more than once in the cluster.
2. If yes: sort the cluster by RT and find the **largest RT gap** between any
   two consecutive peaks.
3. Split the cluster at that gap into two sub-clusters (left and right).
4. Repeat steps 1–3 on each sub-cluster recursively.
5. Terminate when every sub-cluster has unique sample membership, or when a
   sub-cluster contains only one peak (cannot be split further).

**Why split at the largest gap?**
The underlying assumption of RT clustering is that peaks close in RT belong to
the same compound. Splitting at the widest gap is therefore the most conservative
choice: it keeps the most similar peaks together and separates the most
dissimilar ones. It is the natural inverse of the greedy clustering step.

**Outcome:**
After conflict resolution, every cluster has at most one peak per sample.
Every peak is recorded in `feature_peak_log.csv` with `selected = True` — no
data is discarded. Previously, same-sample duplicates would have their
lower-area peak silently dropped; now they each become a distinct feature,
preserving the signal from both compounds.

**Printed output:**
When splits occur, `data_import.py` reports the number of clusters that were
split, e.g.:
```
same-sample conflict splits : 12 cluster(s) split (peaks from the same sample separated into distinct features)
```
A high count suggests a dense chromatogram; consider enabling `USE_MZ = True`
to also separate co-eluters by m/z before the conflict step.

---

### 2.7 Feature ID and metadata

**Feature ID:**
- `USE_MZ = False`: `"RT_{mean_rt:.4f}"` — e.g. `"RT_5.0012"`
- `USE_MZ = True`:  `"{mean_mz:.4f}_{mean_rt:.4f}"` — e.g. `"312.2678_5.0012"`

The mean RT and mean m/z are computed over **all peaks** in the cluster
(including within-sample duplicates), not just the selected ones. This gives a
more robust centroid across samples.

**`feature_metadata.csv` columns:**

| Column | Description |
|--------|-------------|
| `feature_id` | Auto-generated identifier |
| `compound_name` | Name from the most-abundant peak in the cluster |
| `mean_rt` | Mean aligned RT across all cluster peaks (min) |
| `mean_mz` | Mean reference m/z across all cluster peaks (Da) |
| `rt_min`, `rt_max` | RT spread within the cluster |
| `rt_std` | Population SD of RT values within cluster |
| `mz_min`, `mz_max` | m/z spread within the cluster |
| `mz_std` | Population SD of m/z values within cluster |
| `n_samples_detected` | Number of samples where area > 0 |
| `n_samples_total` | Total number of samples |
| `n_contributing_peaks` | Total raw peaks assigned to this cluster |

**`feature_name_map.csv`:** a compact two-column file (`feature_id`,
`compound_name`) for use by plotting modules when `FEATURE_LABEL = "name"`.

**`feature_peak_log.csv` columns:**

| Column | Description |
|--------|-------------|
| `feature_id` | Cluster the peak belongs to |
| `sample` | Source injection file stem |
| `ref_mz` | Reference m/z reported by TraceFinder |
| `rt_raw` | Retention time before alignment |
| `rt_aligned` | Retention time after alignment |
| `rt_shift` | Median shift applied to this sample |
| `area` | Chromatographic peak area (or height) |
| `selected` | `True` if this peak's area entered the matrix |

---

### 2.8 Blank reference table

For each sample feature, the blank reference table records the blank signal that
will be used in Step 2.

**Algorithm (RT-only matching, always performed):**
1. Blank peaks are pooled the same way as sample peaks (§2.3), but alignment is
   applied to the blank files independently with their own reference.
2. For each sample feature, all blank peaks within `±RT_MARGIN` of the
   feature's mean RT are found.
3. The blank peak with the **highest area** within this window is recorded.

Three values are stored per feature:
- `max_blank_area`: the highest blank area within the RT window (0 if no blank
  peak found within the window).
- `blank_rt`: RT of the matched blank peak.
- `blank_mz`: m/z of the matched blank peak.

The m/z of the best-matched blank peak is stored here as audit information
regardless of whether `BLANK_USE_MZ` is enabled — it is used by Step 2's
optional m/z gate.

**Assumption:** Blanks (solvent injections, matrix blanks, or extraction blanks)
represent background signal that is chemically present in the absence of analyte.
A feature that exceeds blank signal by the fold-change threshold is assumed to
be a genuine biological signal.

---

### 2.9 Peak matrix construction

The raw peak matrix is a `features × samples` DataFrame. Rows are feature IDs;
columns are sample names in the order defined by group membership and
alphabetical sort within each group.

- **Detected:** area of the selected peak for that sample.
- **Not detected:** `0` (not `NaN`). This representation is deliberate:
  `log2(0 + 1) = 0`, so undetected features map cleanly to the transformed
  space with no special-case handling needed downstream.

---

## 3. Step 2 — Blank correction

**Script:** `blank_correction.py`
**Input:** `peak_matrix_raw.csv`, `blank_per_feature.csv`, `blank_features.csv`
(fallback), `feature_metadata.csv`, `sample_groups.csv`
**Output:** `peak_matrix_blank_corrected.csv`, `blank_correction_audit.csv`,
`features_removed_blank.csv` (backward-compat summary)

Blank correction removes features — or individual sample cells — that are
plausibly explained by background matrix signal rather than biological analyte.
The background signal is characterised by the blank injections in the sequence.

Two orthogonal dimensions are configurable: how the **sample signal** is
represented (`BLANK_SAMPLE_MODE`) and which **blank signal** is used as the
denominator (`BLANK_REFERENCE_MODE`).

---

### 3.1 Optional m/z gate

**Config:** `BLANK_USE_MZ`, `BLANK_MZ_TOLERANCE`

The blank table from Step 1 performs RT-only matching: any blank peak within
`±RT_MARGIN` of a sample feature is considered a potential match. This creates
a risk of false removal: if a *different compound* (different m/z, different
chemistry) elutes in the blank at the same retention time as a genuine sample
metabolite, the RT-only match would flag the metabolite for removal even though
it is not truly present in the blank.

**When `BLANK_USE_MZ = True`:**

The gate is applied **per blank file** to every row of `blank_per_feature.csv`:
1. Retrieve `mean_mz` of the sample feature from `feature_metadata.csv`.
2. Retrieve `blank_mz` of the RT-matched blank peak.
3. If `|mean_mz − blank_mz| > BLANK_MZ_TOLERANCE`:
   - The blank peak is rejected as a valid match.
   - That blank file's area for this feature is set to `0`.
   - The feature is treated as "not found in that blank file".

If a feature's blank area is zeroed in **all** blank files (every blank fails
the m/z gate), the feature has no valid blank reference and is **always
retained**.  Full gate details — including `mz_delta` and `mz_gate_rejected`
— are recorded in `blank_correction_audit.csv` for every comparison.

**Assumption (m/z gate):** Two compounds are chemically distinct if their
reference m/z values differ by more than `BLANK_MZ_TOLERANCE`. For
high-resolution Orbitrap data, 0.005 Da is a safe threshold that excludes
different elemental compositions while tolerating mass measurement error.

---

### 3.2 Sample-side mode (BLANK_SAMPLE_MODE)

Controls how the sample signal is aggregated for the fold-change comparison.

| Mode | Behaviour | Removal granularity |
|------|-----------|---------------------|
| `"mean"` (default) | Arithmetic mean across **all** samples | Whole feature row removed |
| `"per_sample"` | Each sample compared individually | Failing **sample cells** are zeroed; row dropped only if all cells fail |
| `"per_group"` | Mean per `SAMPLE_GROUPS` group | Failing **group's cells** are zeroed; row dropped only if all groups fail |

**Zero-area skip rule:** If a sample's area is already `0` (feature not
detected in that injection), the blank comparison for that cell is skipped — an
absent signal has no fold change and cannot be further removed.

---

### 3.3 Blank reference mode (BLANK_REFERENCE_MODE)

Controls which blank signal is used as the fold-change denominator.

| Mode | Denominator | Notes |
|------|-------------|-------|
| `"max"` (default) | Highest blank area across all blank files | Most conservative aggregate; matches previous pipeline behaviour |
| `"mean"` | Mean blank area across all blank files | m/z gate applied per file before averaging; a rejected blank contributes 0 to the mean |
| `"each"` | Each blank file compared independently | A sample unit **fails if it fails against any single blank file** (most conservative individual check) |

**Decision rule (per comparison):**
- No valid blank area (blank absent or all m/z-rejected): **always retained**.
- `fold_change ≥ FOLD_CHANGE_THRESHOLD` → **retained**.
- `fold_change < FOLD_CHANGE_THRESHOLD` → **removed** (mean mode) or
  **zeroed** (per_sample / per_group mode).

**Rationale:** A fold-change of 3 means the sample signal is on average 3×
higher than the background. Common practice in untargeted GC-MS metabolomics
uses thresholds of 3 (lenient) to 10 (strict). The choice depends on analyte
class, extraction method, and background complexity.

---

### 3.4 Audit log (blank_correction_audit.csv)

Every comparison — one row per (feature, sample_unit, blank_reference) — is
written to `blank_correction_audit.csv`. This log enables full traceability:

| Column | Description |
|--------|-------------|
| `feature_id` | Feature identifier |
| `mean_rt` | Feature mean RT |
| `mean_mz` | Feature mean m/z |
| `sample_unit` | Sample / group / `"ALL_MEAN"` that was compared |
| `unit_type` | `"sample"` \| `"group"` \| `"all_mean"` |
| `group` | Group name (or `"-"` for all_mean) |
| `sample_area` | Signal area for this comparison unit |
| `blank_name` | Blank file stem, `"BLANK_MAX"`, or `"BLANK_MEAN"` |
| `blank_area` | Blank area after m/z gate |
| `blank_rt` | RT of matched blank peak (NaN if no match) |
| `blank_mz` | m/z of matched blank peak (NaN if no match) |
| `mz_delta` | `|feature_mean_mz − blank_mz|` |
| `mz_gate_rejected` | `True` if blank rejected by m/z gate |
| `fold_change` | `sample_area / blank_area` (NaN = no blank signal) |
| `comparison_failed` | `True` if this specific comparison was below threshold |
| `fold_threshold` | `FOLD_CHANGE_THRESHOLD` applied |
| `decision` | `"kept"` \| `"zeroed"` \| `"removed"` |

`features_removed_blank.csv` is also written for backward compatibility with
tools that read the old pipeline output format.

---

## 4. Step 2b — Prevalence histogram

**Script:** `prevalence_histogram.py`
**Input:** `peak_matrix_blank_corrected.csv`
**Output:** `output/plots/prevalence_histogram.png`, `output/prevalence_summary.csv`

This diagnostic step runs immediately after blank correction and before
normalization. It characterises *how consistently* each feature is detected
across all samples, giving the user a quantitative basis for choosing the
`MIN_PREVALENCE_PCA` and `MIN_PREVALENCE_HCA` thresholds in `config.py`.

---

### 4.1 Prevalence calculation

For each feature (row in the blank-corrected matrix), the script counts the
number of samples in which the feature area is greater than zero:

```
detected_count_i  = Σ_j  [area_ij > 0]
prevalence_i      = detected_count_i / N_samples
```

`prevalence_i` is a value in `{0/N, 1/N, 2/N, …, N/N}` — a discrete set of
`N+1` possible values, where `N` is the total number of sample columns in the
matrix (blanks excluded). A feature with `prevalence = 1.0` is detected in
every sample; `prevalence = 0.0` means the feature was zeroed by blank
correction in all samples and is effectively absent.

---

### 4.2 Bin alignment and x-axis scaling

Because prevalence can only take values at integer multiples of `1/N`, a
standard fixed-width histogram would misalign bars with the true data points
when `N` is not a divisor of 10 or 100. The script instead places bin edges at:

```
edges_k = (k − 0.5) / N    for k = 0, 1, …, N+1
```

This centres each bar exactly on the corresponding `k/N` value. The x-axis
automatically adapts to any number of samples:

- **N ≤ 20:** every tick is labelled as `k/N (xx%)`.
- **N > 20:** ticks are thinned to approximately 10 evenly-spaced labels to
  prevent overplotting, always including the endpoints `0/N` and `N/N`.

Each non-empty bar is annotated with its feature count directly above.

---

### 4.3 Threshold reference lines

Vertical lines are drawn at the prevalence thresholds configured in
`config.py`, shifted to the left edge of the first included bin so the line
marks exactly where the filter cuts:

| Line style | Parameter | Colour |
|------------|-----------|--------|
| Dashed | `MIN_PREVALENCE_PCA` | Red |
| Dotted | `MIN_PREVALENCE_HCA` | Orange |
| Dash-dot | `MIN_PREVALENCE_VOLCANO` | Green |

Lines with a threshold of `0.0` are omitted (no filter applied). This lets
the researcher see at a glance how many features fall below each threshold
and adjust accordingly before committing to normalization.

---

### 4.4 Outputs

**`prevalence_histogram.png`** — bar chart of `# features` (y-axis) vs
`prevalence` (x-axis), one bar per unique `k/N` detection level.

**`prevalence_summary.csv`** — tabular version with columns:

| Column | Description |
|--------|-------------|
| `n_samples_detected` | Integer k — number of samples in which the feature was detected |
| `prevalence` | `k / N_samples` |
| `n_features` | Number of features at exactly this prevalence level |

---

## 5. Step 3 — Normalization, transformation, and scaling

**Script:** `normalization.py`
**Input:** `peak_matrix_blank_corrected.csv`, `feature_metadata.csv`
**Output:** `peak_matrix_processed.csv` (for HCA and volcano),
`peak_matrix_processed_pca.csv` (for PCA)

This step transforms the raw feature areas into a form suitable for multivariate
statistical analysis. Three sequential operations are applied:
normalization → log transformation → scaling. Two output matrices are produced
with different pre-processing filters.

---

### 4.1 Prevalence filter

**Config:** `MIN_PREVALENCE_PCA`, `MIN_PREVALENCE_HCA`, `MIN_PREVALENCE_VOLCANO`

Features that are only detected in a small fraction of samples often represent
noise, carryover, or rare biological events that destabilise multivariate models.
The prevalence filter removes them before the normalization chain.

```
detected_fraction = (n_samples where area > 0) / (total n_samples)
```

A feature is retained if `detected_fraction ≥ MIN_PREVALENCE`.

**Different thresholds for different analyses:**

| Analysis | Recommended default | Rationale |
|----------|--------------------|-----------|
| PCA | 0.5 (≥50% of samples) | Noise features with many zeros destabilise principal components, inflating variance estimates for uninformative dimensions. |
| HCA | 0.0 (disabled) | Sparse features create uninformative columns in the heatmap but do not break the clustering algorithm. Enable if the heatmap is dominated by near-zero columns. |
| Volcano | 0.0 (disabled) | A compound present in 100% of group A but 0% of group B is the most biologically interesting finding — removing it by overall prevalence would be methodologically wrong. |

When removals occur, an audit file is written (`features_removed_prevalence_pca.csv`
or `features_removed_prevalence_hca.csv`) listing the feature_id, number of
detected samples, total samples, and detected fraction.

**Note:** Prevalence filtering is not applied to the volcano analysis within
`normalization.py` (since volcano.py reads the blank-corrected matrix directly
and applies its own optional filter).

---

### 4.2 Exclusion list (PCA matrix only)

**Config:** `EXCLUSION_LIST`, `EXCLUSION_RT_MARGIN`

Some features are already chemically identified and known to dominate the sample
matrix (e.g. a highly abundant internal reference compound or a known artifact).
Including them in PCA can cause the first principal components to describe only
the variance of that one compound, obscuring the biological structure of the
remaining metabolites.

The exclusion list specifies retention times (in minutes) of such compounds.
Any feature whose `mean_rt` falls within `±EXCLUSION_RT_MARGIN` of a listed RT
is excluded **only from the PCA matrix**. It remains in the HCA and volcano
matrices so that it is still statistically testable and visible in those analyses.

Removals are logged to `features_removed_exclusion.csv` with the matched
exclusion RT and the RT deviation.

---

### 4.3 Sum normalization

**Config:** `NORMALIZATION = "sum"`

```
x_normalized(sample j) = x_raw(sample j) / sum_j  ×  median(sum_1, …, sum_n)
```

Where `sum_j = Σ_features(area_ij)` is the total signal in sample `j`.

**Step by step:**
1. Compute the total signal per sample (column sum of the features × samples matrix).
2. Divide each sample by its total signal → all samples now have the same sum
   of 1 (relative composition).
3. Multiply by the median column sum → restores the data to the original
   magnitude scale centred on a typical sample. This prevents downstream
   log-transformation from operating on values that are all near 0 or 1.

**Rationale:** Variation in total signal across samples is often technical
(injection volume, extraction yield, instrument drift) rather than biological.
Sum normalization removes this global scaling factor. It is the standard
first-choice normalization for untargeted GC-MS data.

**Assumption:** The total metabolic output (total ion current) does not differ
biologically between groups. If a treatment genuinely increases or decreases
overall metabolism, sum normalization will suppress this information. In that
case, an external reference (e.g. an internal standard or a dilution series)
should be used instead.

---

### 4.4 Median normalization (alternative)

**Config:** `NORMALIZATION = "median"`

```
x_normalized(sample j) = x_raw(sample j) / median_j  ×  global_median
```

Where `median_j` is the median non-zero feature intensity in sample `j`, and
`global_median` is the median of all sample medians.

This is more robust than sum normalization when many features are zero (zeros
drag the mean downward but do not affect the median). Useful when most features
are absent in most samples (sparse matrices).

---

### 4.5 Log transformation

**Config:** `LOG_BASE`

```
x_log = log2(x + 1)
```

The pseudo-count of +1 is added before taking the logarithm so that undetected
features (area = 0) map cleanly to 0 rather than −∞. This is the standard
approach in metabolomics and genomics.

**Why log2?**
- Mass spectrometry signal intensities span several orders of magnitude
  (typically 3–6 decades). Without transformation, high-abundance features
  completely dominate Euclidean-distance-based analyses.
- Log-transformation compresses dynamic range and stabilises variance across
  features, moving the distribution toward approximate normality.
- Log2 is conventional in metabolomics: a difference of 1 unit = a 2-fold
  difference in abundance, which is directly interpretable.
- A difference of 1 in log2 space equals a 2-fold change.

**Assumption:** Signal intensities are approximately log-normally distributed —
i.e. the relative rather than absolute differences between samples are
biologically relevant. This is a well-validated assumption for GC-MS peak areas.

---

### 4.6 Pareto scaling (recommended)

**Config:** `SCALING = "pareto"`

Applied per feature (row-wise across all samples):

```
x_scaled = (x - mean_feature) / sqrt(std_feature)
```

Where `mean_feature` and `std_feature` (sample SD, ddof=1) are computed across
all samples for that feature.

**What it does:**
- Mean-centring (`x − mean`): moves the data to be centred on 0. Required for
  PCA and HCA to function correctly (otherwise the first PC captures only
  the mean, not variance structure).
- Division by `sqrt(std)`: partially reduces the dominance of high-variance
  (typically high-abundance) features, but less aggressively than auto-scaling.

**Why Pareto over auto-scaling?**
- Auto-scaling (division by `std`) gives every feature exactly unit variance.
  This is mathematically elegant but treats a feature with a tiny, noisy signal
  the same as a feature with a large, reproducible signal. The result is that
  noise features can dominate PCA components.
- Pareto scaling retains more biological variance structure by only partially
  down-weighting high-abundance features.
- For untargeted GC-MS, Pareto is the consensus recommendation (van den Berg
  et al. 2006, BMC Genomics).

**Zero-variance features:** Features that are identical across all samples after
log-transformation (std = 0) cannot be scaled and carry no discriminating
information. They are automatically dropped and logged to
`features_removed_zero_variance.csv`.

---

### 4.7 Auto scaling (alternative)

**Config:** `SCALING = "auto"`

```
x_scaled = (x - mean_feature) / std_feature
```

Gives every feature exactly unit variance after scaling. Use when all features
should contribute equally to the analysis regardless of abundance, for example
in targeted analyses where all analytes are of equal interest.

---

### 4.8 Two output matrices

Both matrices are transposed to `samples × features` (rows = samples, columns
= features) — the standard convention for scikit-learn and most ML/stats libraries.

| Matrix | Filters applied | Used by |
|--------|----------------|---------|
| `peak_matrix_processed.csv` | Prevalence (HCA threshold), then normalize + log + scale | HCA (Step 5), Volcano (Step 6) via re-normalization |
| `peak_matrix_processed_pca.csv` | Prevalence (PCA threshold) + exclusion list, then normalize + log + scale | PCA (Step 4) |

When `MIN_PREVALENCE_HCA = 0` and `EXCLUSION_LIST = []` and the two prevalence
thresholds are identical, the two matrices are identical and the PCA matrix is
a reference to the same object (not duplicated in memory).

---

## 6. Step 4 — Principal Component Analysis (PCA)

**Script:** `pca.py`
**Input:** `peak_matrix_processed_pca.csv`, `sample_groups.csv`
**Output:** Scores CSV, loadings CSV, variance CSV, PNG plots, optional 3D HTML

### 5.1 Algorithm

PCA is computed using `sklearn.decomposition.PCA` with `n_components =
min(N_COMPONENTS, n_samples, n_features)` (truncated to the rank of the
input matrix if necessary).

Scikit-learn's PCA uses full Singular Value Decomposition (SVD):

```
X = U Σ Vᵀ
```

Where `X` is the mean-centred data matrix (samples × features), `U` is the
matrix of left singular vectors (sample scores, unnormalised), `Σ` is the
diagonal matrix of singular values, and `Vᵀ` holds the right singular vectors
(principal component directions, i.e. loadings).

**Scores** (`pca_scores.csv`):
`scores = X × V` — the projection of each sample onto each PC axis.
Each row is a sample; each column is a PC. Units are arbitrary (scaled by
the data's variance structure).

**Loadings** (`pca_loadings.csv`):
`loadings = V` — the weight of each original feature in each PC. Stored as
`pca.components_.T` (shape: features × n_components). A large positive loading
on PC1 means the feature contributes strongly to positive PC1 scores.

**Explained variance** (`pca_variance.csv`):
- `explained_variance`: eigenvalue (variance along each PC in the original
  feature space).
- `explained_variance_ratio`: fraction of total dataset variance explained.
- `cumulative_explained_variance`: running sum of ratios.

**Why PCA on normalised, log2-transformed, Pareto-scaled data?**
PCA is a linear method that maximises explained variance. It is sensitive to the
scale of features: without normalization and scaling, one high-abundance
compound would dominate all PCs. The preprocessing pipeline (§4.3–4.6) brings
all features to a comparable scale before PCA is applied.

**Assumption:** The main sources of biological variation in the data are linear
combinations of the measured features. Non-linear structures (e.g. batch
effects, non-monotone dose-responses) may require kernel PCA or t-SNE instead.

---

### 5.2 Scores plot and confidence ellipses

**Config:** `PCA_PLOT_X`, `PCA_PLOT_Y`, `PCA_ELLIPSE`

The scores plot is a scatter plot of PC_x vs PC_y for all samples, coloured by
group. Sample names are annotated as text offsets.

**Confidence ellipses (`PCA_ELLIPSE = True`):**

For each group with ≥2 samples, a 95% confidence ellipse is drawn using
eigendecomposition of the 2×2 covariance matrix of the group's PC coordinates:

```
C = [[cov(PCx, PCx), cov(PCx, PCy)],
     [cov(PCy, PCx), cov(PCy, PCy)]]
```

The ellipse parameters are derived from the eigenvectors and eigenvalues of C:
- **Centre:** `(mean(PC_x), mean(PC_y))` for the group.
- **Semi-axes:** `n_std × sqrt(eigenvalue_i)` for each eigenvector direction,
  with `n_std = 2.0`.
- **Rotation angle:** `arctan2` of the first eigenvector.

`n_std = 2.0` produces an ellipse that contains approximately 86% of a
bivariate normal distribution (not 95%; the 95% level for bivariate normal
requires `n_std ≈ 2.45`). The value is kept at 2.0 as a practical visual aid
rather than a strict statistical boundary.

**Limitation:** The ellipses are parametric and assume bivariate normality of
the group's sample scores. With small sample sizes (n < 5), they are
unreliable guides and should be treated as visual orientation, not inference.

---

### 5.3 Loadings plot and bar chart

**Loadings scatter (`pca_loadings.png`):**
Each feature is plotted as a point at `(loading_PC_x, loading_PC_y)`. The top
`PCA_TOP_LOADINGS` features by Euclidean distance from the origin are labelled
and highlighted in red.

```
distance_i = sqrt(loading_PC_x_i² + loading_PC_y_i²)
```

Features far from the origin drive the separation of samples along those PCs.

**Loading bar chart (`pca_loadings_bar.png`):**
The same top-N selection is used. Features are sorted by their PC_x loading
(ascending) so that the direction of effect (positive vs negative PC_x) is
immediately visible. Each feature gets two bars — one for PC_x and one for PC_y
— shown simultaneously, with value labels at bar ends.

**Feature labels in plots:**
When `FEATURE_LABEL = "name"`, compound names from `feature_name_map.csv` are
used as axis labels and annotations in place of feature IDs.

---

## 7. Step 5 — Hierarchical Cluster Analysis (HCA)

**Script:** `hca.py`
**Input:** `peak_matrix_processed.csv` (full feature set), `sample_groups.csv`
**Output:** `hca_heatmap.png`, `hca_row_order.csv` (sample order), `hca_col_order.csv` (feature order)

### 6.1 Algorithm and linkage

HCA builds a hierarchy of clusters by iteratively merging the two most similar
clusters. The result is displayed as a dendrogram; the heatmap columns and rows
are reordered to follow the dendrogram leaf order, visually grouping similar
samples and similar features together.

**Distance metric (`HCA_METRIC`):**

| Option | Formula | Use case |
|--------|---------|----------|
| `euclidean` | `sqrt(Σ(xi − yi)²)` | Default; required for Ward linkage |
| `correlation` | `1 − Pearson_r(xi, yi)` | Focuses on pattern shape, not magnitude |

**Linkage method (`HCA_LINKAGE`):**

| Method | Merge criterion | Properties |
|--------|----------------|------------|
| `ward` (default) | Minimise increase in total within-cluster sum of squares | Tends to produce compact, equally-sized clusters; preferred for metabolomics |
| `average` | Average pairwise distance between all inter-cluster pairs | Robust, less biased toward cluster size |
| `complete` | Maximum pairwise distance (furthest-neighbour) | Forms tight, spherical clusters |
| `single` | Minimum pairwise distance (nearest-neighbour) | Susceptible to chaining effect |

**Ward linkage in detail:**
At each step, the algorithm merges the two clusters `A` and `B` that minimise
the **Ward criterion** (increase in error sum of squares):

```
ΔW(A, B) = (|A| × |B|) / (|A| + |B|) × ||mean_A − mean_B||²
```

This is equivalent to minimising the within-cluster variance after merging.
Ward linkage can only be used with Euclidean distance (it requires the
Euclidean geometry to define "within-cluster variance").

**Assumption:** Samples that cluster together are more similar in their overall
metabolic profile than samples in different clusters. The clustering identifies
co-varying groups of features and separates biologically distinct sample classes
without requiring pre-defined labels.

---

### 6.2 Heatmap construction

The heatmap is produced by `seaborn.clustermap`, which computes hierarchical
clustering on both rows (samples) and columns (features) simultaneously.

**Colour mapping:**
- Colormap `"vlag"` (default): a diverging blue-white-red palette.
- `center = 0`: the colourmap is centred at 0, which is correct for
  Pareto-scaled (mean-centred) data. Blue = below mean; red = above mean.
- A feature consistently red in one group and blue in another is a potential
  biomarker.

**Row colour annotation:** A colour strip on the left of the heatmap encodes
sample group membership, using the same palette as the PCA scores plot.

**Feature labels:**
Column (feature) labels are shown when `n_features ≤ HCA_MAX_FEATURE_LABELS`.
When `FEATURE_LABEL = "name"`, compound names are used as column labels.

**Dendrogram orders exported:**
`hca_row_order.csv` and `hca_col_order.csv` record the integer leaf order from
the row and column dendrograms respectively. This order is useful for:
- Extracting co-varying metabolite clusters for downstream Spearman correlation.
- Reporting the exact cluster order in supplementary tables.
- Reproducing the figure ordering in other tools.

---

## 8. Step 6 — Volcano plot

**Script:** `volcano.py`
**Input:** `peak_matrix_blank_corrected.csv` (raw blank-corrected areas),
`sample_groups.csv`
**Output:** Per comparison: `volcano_{A}_vs_{B}.png`, `volcano_{A}_vs_{B}.csv`

The volcano plot is a per-feature scatter plot of `log2 fold change` (x-axis)
vs `−log10(adjusted p-value)` (y-axis). Features with large fold changes and
high statistical significance appear in the upper corners and are the primary
candidates for biological relevance.

---

### 7.1 Normalization for volcano (no Pareto scaling)

The blank-corrected raw matrix (not the Pareto-scaled matrix from Step 3) is
used as input. Within `volcano.py`, normalization and log transformation are
applied inline using the same formulae as Step 3 (§4.3, §4.5), but **Pareto
scaling is intentionally omitted**.

**Why no Pareto scaling for volcano?**
Pareto scaling mean-centres each feature and divides by the square root of its
standard deviation. This changes the *absolute values* and *relative magnitudes*
between features. Fold change (`mean_A / mean_B`) computed from Pareto-scaled
values would not reflect the true biological fold difference — it would reflect
a scaled, mean-centred ratio that has no direct biological interpretation.
Using unscaled log-normalised values preserves the quantitative meaning of
fold change.

---

### 7.2 Log2 fold change

For each feature:

```
log2FC = mean(log2(area + 1) in group A)  −  mean(log2(area + 1) in group B)
```

This is the **difference of log2 means**, which equals the **log2 of the ratio
of geometric means**:

```
log2FC = log2(geometric_mean_A / geometric_mean_B)
```

**Interpretation:**
- `log2FC = 1.0` → group A is 2× higher than group B (2-fold increase).
- `log2FC = −1.0` → group A is 2× lower (group B is higher).
- `log2FC = 0` → no difference.
- The sign convention is always `A / B`: positive means higher in group A.

**Sign convention note:** When `VOLCANO_COMPARISONS = "all"`, all pairwise
comparisons are run in alphabetical order. If you specify comparisons manually
as `[("treatment", "control")]`, the sign will indicate treatment vs control.

**Effect of zeros:**
`log2(0 + 1) = 0`. When a feature is absent from all samples in one group
(area = 0 for all samples), its mean log2 for that group is 0. The fold change
will then equal the mean log2 of the other group, correctly reflecting a
group-specific feature. This is why `MIN_PREVALENCE_VOLCANO = 0.0` is the
default: these group-specific features often have the largest fold changes and
the strongest biological significance.

---

### 7.3 Mann-Whitney U test

For each feature, the Mann-Whitney U test (also called the Wilcoxon rank-sum
test) assesses whether the distribution of log-normalised values in group A is
stochastically greater than in group B.

**Why Mann-Whitney?**
- Non-parametric: makes no assumption about the distribution of feature
  intensities. With typical metabolomics sample sizes (n = 3–10 per group),
  normality cannot be verified and a t-test would be unreliable.
- Robust to outliers: based on ranks rather than raw values.
- Handles ties: the scipy implementation uses a continuity correction.
- Two-sided alternative: tests for any difference in location (not just in
  one direction).

**Degenerate case handling:**
If all values in one group are identical (e.g. all zero), scipy raises a
`ValueError` because ranks cannot be assigned. This is caught and the p-value
is set to `1.0` (no evidence of difference).

**Limitation:**
Mann-Whitney tests whether the distributions are *stochastically ordered* but
does not directly test the mean difference. With very small groups (n = 2–3),
the minimum achievable p-value is bounded away from 0 regardless of the effect
size (e.g. with n=3 in each group, the minimum two-sided p-value is ~0.1).

---

### 7.4 Benjamini-Hochberg FDR correction

**Config:** `VOLCANO_P_THRESHOLD`

When testing many features simultaneously (commonly hundreds to thousands in
untargeted metabolomics), the probability of at least one false positive grows
rapidly. Controlling the per-test Type-I error rate at α = 0.05 for 500
features would yield ~25 expected false positives.

The Benjamini-Hochberg (BH) procedure controls the **False Discovery Rate (FDR)** —
the expected proportion of false discoveries among all rejected hypotheses.

**Algorithm:**

1. Rank the `m` raw p-values in ascending order: `p_(1) ≤ p_(2) ≤ … ≤ p_(m)`.
2. Compute adjusted p-values:
   ```
   p_adj(i) = p_raw(i) × m / rank(i)
   ```
3. Clip to [0, 1].
4. Enforce monotonicity (backward pass): for `i` from `m−1` to `1`:
   ```
   p_adj_sorted(i) = min(p_adj_sorted(i), p_adj_sorted(i+1))
   ```
   This ensures that a feature with a lower raw p-value cannot have a higher
   adjusted p-value than a feature with a higher raw p-value.

**Interpretation:**
A feature with `adj_pvalue = 0.03` at a FDR threshold of 0.05 means: among
all features called significant at this threshold, we expect ≤5% to be false
positives. It does **not** mean the probability of that specific feature being
a false positive is 3%.

**Comparison to Bonferroni:**
Bonferroni controls FWER (family-wise error rate, probability of *any* false
positive) by multiplying p-values by `m`. BH is less stringent and more
statistically powerful, accepting a small proportion of false discoveries in
exchange for a lower false-negative rate. BH is the standard choice in
metabolomics.

---

### 7.5 Classification and thresholds

**Config:** `VOLCANO_FC_THRESHOLD`, `VOLCANO_P_THRESHOLD`

Each feature is assigned one of three labels:

| Label | Condition |
|-------|-----------|
| `up_A` | `log2FC > FC_THRESHOLD`  AND  `adj_pvalue < P_THRESHOLD` |
| `up_B` | `log2FC < −FC_THRESHOLD`  AND  `adj_pvalue < P_THRESHOLD` |
| `n.s.` | All other features |

**Threshold interpretation:**
- `VOLCANO_FC_THRESHOLD = 1.0` → 2-fold change minimum (log2(2) = 1.0). A
  common choice for biological significance in metabolomics.
- `VOLCANO_P_THRESHOLD = 0.05` → 5% FDR. Standard threshold; tighten to 0.01
  for greater confidence when sample size permits.

The `VOLCANO_TOP_LABELS` most significant (lowest adj_pvalue) features that
exceed both thresholds are labelled in the plot with their feature ID or compound
name (when `FEATURE_LABEL = "name"`).

The full results table (`volcano_{A}_vs_{B}.csv`) contains all features with
their `log2FC`, `pvalue`, `adj_pvalue`, and `direction`, regardless of
significance. This allows the user to re-apply different thresholds without
re-running the analysis.

---

## 9. Report — Blank contaminants

**Script:** `blank_contaminants_report.py`
**Input:** `blank_correction_audit.csv`, `feature_name_map.csv`
**Output:** `blank_contaminants_report.csv`

This read-only reporting step filters the full audit log from Step 2 to show
only comparisons where a feature was **removed** or a sample cell was
**zeroed**. Each output row corresponds to one (feature, sample_unit, blank)
comparison, so every removal decision can be traced back to the exact blank
file and fold change that triggered it.

**Output columns:**

| Column | Description |
|--------|-------------|
| `feature_id` | Feature identifier |
| `compound_name` | Compound name from `feature_name_map.csv` (empty if unknown) |
| `mean_rt` | Feature mean aligned RT (4 dp) |
| `mean_mz` | Feature mean reference m/z (4 dp) |
| `sample_unit` | Sample name / group name / `"ALL_MEAN"` |
| `unit_type` | `"sample"` \| `"group"` \| `"all_mean"` |
| `group` | Group the sample_unit belongs to |
| `sample_area` | Signal area for this comparison unit |
| `blank_name` | Blank file stem, `"BLANK_MAX"`, or `"BLANK_MEAN"` |
| `blank_area` | Blank area used (after m/z gate) |
| `blank_rt` | RT of matched blank peak |
| `blank_mz` | m/z of matched blank peak |
| `mz_delta` | `|feature_mean_mz − blank_mz|` |
| `mz_gate_rejected` | `True` if blank was rejected by the m/z gate |
| `fold_change` | `sample_area / blank_area` |
| `comparison_failed` | `True` if this comparison individually was below threshold |
| `fold_threshold` | `FOLD_CHANGE_THRESHOLD` applied |
| `decision` | `"removed"` (feature fully dropped) or `"zeroed"` (cell set to 0) |

**How to read the report with `BLANK_REFERENCE_MODE = "each"`:**
In this mode, one row appears per (feature, sample_unit, blank_file). A
feature/cell is removed when at least one blank triggers failure.
`comparison_failed = True` identifies the specific blank(s) responsible.

**Backward compatibility:**
If `blank_correction_audit.csv` is absent (old pipeline run), the script falls
back to reading `features_removed_blank.csv` and `feature_peak_log.csv` and
produces the legacy 6-column format.

**Why is this useful?**
The blank-corrected matrix simply records zeros for removed features — there is
no record in the matrix itself of *why* a feature was removed or *which* blank
triggered it. This report restores full provenance:
- Siloxanes (column bleed) and solvent peaks are typical removed features.
- In `per_sample` or `per_group` mode, you can see that a feature was removed
  from one group's cells but retained for another — revealing group-specific
  blank contamination.
- `mz_gate_rejected` rows show where m/z-mismatched blank peaks were
  discarded, protecting genuine sample metabolites from false removal.

---

## 10. Report — Top PCA features

**Script:** `top_features_analysis.py`
**Input:** `pca_loadings.csv`, `peak_matrix_raw.csv`, `feature_metadata.csv`,
raw TraceFinder CSV files in `DATA/`
**Output:** `top_features_analysis.csv`

This read-only reporting step identifies the features with the highest influence
on the PCA result and summarises their peak areas per sample together with
compound names resolved from the raw data.

**Feature selection — loading magnitude:**

For each feature, the Euclidean magnitude of its loading vector is computed:

$$\text{magnitude}_i = \sqrt{PC1_i^2 + PC2_i^2 (+ PC3_i^2)}$$

PC3 is included only when a third component was computed (`N_COMPONENTS = 3`).
The top `N` features (default 10, configurable via `--num-features`) with the
highest magnitude are selected. These are the features that most strongly drive
sample separation in the PCA scores plot.

**Compound name resolution:**
For each top feature, the script searches the raw TraceFinder CSV files for
rows whose `Retention Time` falls within a ±0.1 min window of the feature's
`mean_rt`. The `Component Name` is read from whichever sample file has the
best (closest RT) match among samples where the feature was detected. The
same `skiprows` auto-detection logic used in `data_import.py` is applied so
that TraceFinder files with an extra metadata header row are handled correctly.

**Output columns:**

| Column | Description |
|--------|-------------|
| `feature_id` | Feature identifier |
| `compound_name` | Name from raw data by RT matching ("Unknown (RT x.xxxx)" if not found) |
| `RT` | Feature mean retention time |
| `area` | Raw peak area for this sample |
| `sample` | Sample name |
| `PC1` | Feature loading on PC1 |
| `PC2` | Feature loading on PC2 |
| `PC3` | Feature loading on PC3 (column absent when `N_COMPONENTS < 3`) |

One row is written per detected sample per feature (only samples where area > 0
are included). Rows are sorted by loading magnitude (descending) then area
(descending).

**Standalone usage:**
```bash
python top_features_analysis.py          # top 10 features
python top_features_analysis.py --n 20  # top 20 features
```

---

## 11. Feature labelling (compound names)

**Config:** `COMPOUND_NAME_COL`, `FEATURE_LABEL`

By default, features are labelled by their auto-generated ID (e.g.
`"RT_5.0012"`), which encodes the mean retention time. When `FEATURE_LABEL = "name"`,
compound names are used in plot annotations, axis tick labels, and bar chart
y-axes.

**How compound names are assigned:**

1. During Step 1, the column specified by `COMPOUND_NAME_COL` (default:
   `"Component Name"`) is read from each TraceFinder CSV row.
2. Within each feature cluster, the peak with the **highest area across all
   samples** is identified. Its compound name is assigned as the feature's
   name.
3. The name is stored in `feature_metadata.csv` (column `compound_name`) and
   in the convenience file `feature_name_map.csv`.

**Why the highest-area peak?**
The most abundant detection is the most reliable. TraceFinder assigns compound
names based on library matching; the highest-area peak typically has the best
signal-to-noise ratio and thus the most confident library match.

**Fallback:** If `compound_name` is empty or NaN (e.g. the TraceFinder CSV had
no name column, or the peak was unidentified), the feature_id is used as the
label. This ensures plots are always fully annotated.

**Output CSVs are not affected:** Regardless of `FEATURE_LABEL`, all output CSV
files use `feature_id` as the index. This ensures that CSVs can always be joined
back to `feature_metadata.csv` and `feature_peak_log.csv` unambiguously.

---

## 12. Assumptions summary

The following is a consolidated list of the key assumptions embedded in the
pipeline design. Violating these assumptions does not always produce wrong
results, but the user should be aware of them when interpreting output.

### Instrument and acquisition

| Assumption | Where used | Consequence if violated |
|-----------|-----------|------------------------|
| A single median RT shift per sample corrects all RT drift | RT alignment (§2.2) | Compounds with RT-dependent drift (non-uniform across the chromatogram) may be split into multiple features |
| `"Reference m/z"` is a stable, compound-specific identifier | RT alignment, m/z clustering, blank m/z gate | Ion switching or calibration errors would misalign peaks |
| Peak area is proportional to analyte concentration over the measured range | All quantitative steps | Saturation or ion suppression breaks this linearity |
| All injections belong to the same analytical sequence (same column, carrier gas, temperature programme) | RT alignment | If comparing across sequences, RT drift may exceed `RT_MARGIN` |

### Feature detection

| Assumption | Where used | Consequence if violated |
|-----------|-----------|------------------------|
| The same compound elutes within `RT_MARGIN` minutes across all samples after alignment | Feature clustering (§2.4) | Features split across samples are counted as different features; false zero entries in the matrix |
| Co-eluting compounds are merged into one feature regardless of their m/z when `USE_MZ = False` | m/z sub-clustering (§2.5) | All co-eluters (including those with different m/z that `USE_MZ = True` would separate) are merged; the area reflects the combined signal |
| The largest RT gap between same-sample peaks marks the boundary between two distinct features | Conflict resolution (§2.6) | If two different compounds from the same sample have very similar RTs, the split may not cleanly separate them; enabling `USE_MZ = True` resolves this by adding m/z as a second criterion |

### Blank correction

| Assumption | Where used | Consequence if violated |
|-----------|-----------|------------------------|
| Blank injections faithfully represent the background matrix signal | Fold-change filter (§3.2) | If blanks were injected under different conditions, background may be over- or under-estimated |
| The maximum blank area is the most conservative (worst-case) background estimate | Fold-change filter | If one blank is an outlier (contaminated injection), genuine features may be wrongly removed |
| Different m/z at the same RT implies a different compound (`BLANK_USE_MZ = True`) | m/z gate (§3.1) | True if mass accuracy is high; less reliable with low-resolution data |

### Normalization and transformation

| Assumption | Where used | Consequence if violated |
|-----------|-----------|------------------------|
| Differences in total signal are technical, not biological | Sum normalization (§4.3) | If total metabolic activity genuinely differs between groups, sum normalization suppresses a real biological signal |
| Feature intensities are approximately log-normally distributed | Log2 transformation (§4.5) | If the distribution is multimodal or highly skewed even after log, parametric methods downstream are unreliable |
| Pareto scaling does not destroy the biological variance structure | Pareto scaling (§4.6) | For very sparse data (many zeros), mean-centring can introduce artefactual structure |

### Multivariate analysis

| Assumption | Where used | Consequence if violated |
|-----------|-----------|------------------------|
| Biological variation is captured in a linear combination of features | PCA (§5.1) | Non-linear metabolic relationships are not captured; use kernel PCA or UMAP instead |
| Samples within a group are biologically homogeneous (no batch effects within group) | PCA ellipses, HCA, volcano | Batch effects appear as within-group separation that mimics biological structure |
| Ward linkage Euclidean distance is an appropriate similarity metric for these data | HCA (§6.1) | For non-Euclidean similarity (e.g. correlation patterns), use "average"/"complete" linkage with "correlation" metric |
| Mann-Whitney U can detect meaningful differences in small groups | Volcano (§7.3) | With n=3 per group, power is low; only large effect sizes will reach significance |
| BH FDR correction is appropriate (tests are not strongly positively correlated) | Volcano (§7.4) | If features are highly correlated (common in metabolomics pathway clusters), BH may be slightly anti-conservative |

---

## 13. Key parameter reference

| Parameter | Default | Step | Description |
|-----------|---------|------|-------------|
| `BLANK_PREFIX` | `"Blank"` | 1 | Filename prefix identifying blank injections |
| `SAMPLE_GROUPS` | — | 1 | Ordered list of `(group_name, prefix)` pairs |
| `RT_MARGIN` | `0.05` | 1 | Half-window (min) for RT feature clustering |
| `USE_MZ` | `False` | 1 | Also cluster by m/z to separate co-eluters |
| `MZ_TOLERANCE` | `0.005` | 1 | Da — m/z proximity for feature-level clustering |
| `ALIGN_RT` | `True` | 1 | Apply median RT-shift correction across samples |
| `MZ_ALIGN_TOLERANCE` | `0.1` | 1 | Da — m/z window for RT alignment matching |
| `VALUE_COL` | `"Area"` | 1 | Column to use: `"Area"` or `"Height"` |
| `COMPOUND_NAME_COL` | `"Component Name"` | 1 | Column holding compound names in TraceFinder CSVs |
| `FOLD_CHANGE_THRESHOLD` | `3.0` | 2 | Min sample/blank area ratio to retain a feature |
| `BLANK_USE_MZ` | `False` | 2 | Also gate blank matching by m/z proximity |
| `BLANK_MZ_TOLERANCE` | `0.005` | 2 | Da — max \|feature_mz − blank_mz\| for blank match |
| `BLANK_SAMPLE_MODE` | `"mean"` | 2 | `"mean"` \| `"per_sample"` \| `"per_group"` — how to aggregate samples for blank comparison |
| `BLANK_REFERENCE_MODE` | `"max"` | 2 | `"max"` \| `"mean"` \| `"each"` — which blank signal to use as the fold-change denominator |
| `MIN_PREVALENCE_PCA` | `0.5` | 3 | Min fraction of samples detecting a feature (for PCA) |
| `MIN_PREVALENCE_HCA` | `0.0` | 3 | Same filter for HCA (0.0 = disabled) |
| `MIN_PREVALENCE_VOLCANO` | `0.0` | 3 | Same filter for volcano (keep at 0.0) |
| `EXCLUSION_LIST` | `[]` | 3 | RT values (min) of known features to exclude from PCA |
| `EXCLUSION_RT_MARGIN` | `0.05` | 3 | ±window for exclusion list RT matching |
| `NORMALIZATION` | `"sum"` | 3 | `"sum"` \| `"median"` \| `"none"` |
| `LOG_BASE` | `2` | 3 | Log base: `2` (log2) or `math.e` (natural log) |
| `SCALING` | `"pareto"` | 3 | `"pareto"` \| `"auto"` \| `"none"` |
| `N_COMPONENTS` | `2` | 4 | Number of PCA components to compute |
| `PCA_PLOT_X` | `1` | 4 | PC for x-axis of 2D plot (1-indexed) |
| `PCA_PLOT_Y` | `2` | 4 | PC for y-axis of 2D plot (1-indexed) |
| `PCA_ELLIPSE` | `True` | 4 | Draw 95% confidence ellipses per group |
| `PCA_TOP_LOADINGS` | `10` | 4 | Features labelled in loadings scatter plot |
| `PCA_BAR_TOP` | `10` | 4 | Features shown in loadings bar chart |
| `HCA_LINKAGE` | `"ward"` | 5 | Hierarchical clustering linkage method |
| `HCA_METRIC` | `"euclidean"` | 5 | Distance metric for HCA |
| `HCA_CMAP` | `"vlag"` | 5 | Colormap for heatmap |
| `HCA_MAX_FEATURE_LABELS` | `50` | 5 | Show feature labels when n_features ≤ this |
| `VOLCANO_COMPARISONS` | `"all"` | 6 | `"all"` or list of `(groupA, groupB)` tuples |
| `VOLCANO_FC_THRESHOLD` | `1.0` | 6 | log2 fold-change cutoff (1.0 = 2-fold) |
| `VOLCANO_P_THRESHOLD` | `0.05` | 6 | BH-adjusted p-value significance threshold |
| `VOLCANO_TOP_LABELS` | `10` | 6 | Most significant features labelled in plot |
| `FEATURE_LABEL` | `"id"` | 4,5,6 | `"id"` = feature_id in plots; `"name"` = compound name |

---

*Generated for pipeline version as of 2026-03-31. All formulae are derived
directly from the Python source code; parameters refer to `config.py`.*
