# TF_NT - Nontargeted GC-MS Metabolomics Pipeline for TraceFinder CSV exports

A modular, configurable data processing pipeline for nontargeted GC-MS metabolomics
data exported from **Thermo TraceFinder** (Orbitrap MS).

## What it does

| Step | Script | Input -> Output |
|------|--------|----------------|
| 1 | `data_import.py` | TraceFinder CSVs -> `peak_matrix_raw.csv` |
| 2 | `blank_correction.py` | raw matrix -> `peak_matrix_blank_corrected.csv` |
| 2b | `prevalence_histogram.py` | `peak_matrix_blank_corrected.csv` -> `prevalence_histogram.png` + `prevalence_summary.csv` |
| 2c | `compound_classification.py` | `feature_metadata.csv` -> `compound_classes.csv` + `feature_metadata_enriched.csv` *(optional; see `RUN_COMPOUND_CLASSIFICATION`)* |
| 2d | `compound_class_plots.py` | `feature_metadata_enriched.csv` -> `class_pie_*.png` *(optional; see `RUN_CLASS_PLOTS`)* |
| 3 | `normalization.py` | blank-corrected -> `peak_matrix_processed.csv` + `peak_matrix_processed_pca.csv` |
| 4 | `pca.py` | `peak_matrix_processed_pca.csv` -> scores/loadings plots & CSVs |
| 5 | `hca.py` | `peak_matrix_processed.csv` -> clustered heatmap + dendrogram order CSVs |
| 6 | `volcano.py` | `peak_matrix_blank_corrected.csv` -> volcano plots & results tables |
| 7 | `top_features_analysis.py` | `pca_loadings.csv` + raw CSVs -> `top_features_analysis.csv` |
| 8 | `blank_contaminants_report.py` | `features_removed_blank.csv` + provenance files -> `blank_contaminants_report.csv` |

All outputs are written to the `output/` folder.

### Note on the exclusion list

`EXCLUSION_LIST` in `config.py` lists retention times of biologically relevant
features to withhold from PCA (so they do not dominate the principal components).
It is applied only when building `peak_matrix_processed_pca.csv` (Step 3 -> Step 4).
HCA and the volcano plot always use the full feature set so that known compounds
remain statistically testable and visible in those analyses.

## Setup

### Requirements

- Python 3.9 or newer

### Installation

```bash
# 1. Create a virtual environment
python3 -m venv .venv

# 2. Activate it
source .venv/bin/activate        # macOS / Linux
.venv\Scripts\activate           # Windows

# 3. Install dependencies
pip install -r requirements.txt
```

## Usage

### 1. Add your data

Place TraceFinder CSV exports in the `DATA/` folder.
See `DATA/README.txt` for expected file naming and column format.

### 2. Configure the pipeline

Open `config.py` and adjust the parameters for your experiment:

```python
# Key settings to review for every new dataset:
BLANK_PREFIX  = "Blank"          # prefix of blank file names
SAMPLE_GROUPS = [                # define your groups (most specific prefix first)
    ("S-R", "S-R"),
    ("S",   "S"),
]
RT_MARGIN             = 0.05     # RT clustering window (minutes)
FOLD_CHANGE_THRESHOLD = 3.0      # blank correction stringency
EXCLUSION_LIST        = []       # RT values (min) of features to exclude from PCA only
NORMALIZATION         = "sum"    # "sum" | "median" | "none"
SCALING               = "pareto" # "pareto" | "auto" | "none"
N_COMPONENTS          = 2        # PCA components (use 3 for interactive 3D plots)
VOLCANO_COMPARISONS   = "all"    # "all" for all pairwise, or [(groupA, groupB), ...]
```

### 3. Run the pipeline

```bash
python pipeline.py
```

Or run individual steps independently (useful when re-running after parameter changes):

```bash
python data_import.py           # Step 1  - re-run if data or RT_MARGIN changes
python blank_correction.py      # Step 2  - re-run if FOLD_CHANGE_THRESHOLD changes
python prevalence_histogram.py      # Step 2b - re-run after blank correction to review detection rates
python compound_classification.py  # Step 2c - re-run to refresh PubChem annotations (uses cache)
python compound_class_plots.py     # Step 2d - re-run after classification or if CLASS_PIE_* settings change
python normalization.py            # Step 3  - re-run if NORMALIZATION, SCALING, or EXCLUSION_LIST changes
python pca.py                # Step 4 - re-run if PCA settings change
python hca.py                # Step 5 - re-run if HCA_LINKAGE or HCA_METRIC changes
python volcano.py            # Step 6 - re-run if VOLCANO_* thresholds change
python top_features_analysis.py          # Step 7 - re-run after PCA; optional --n flag
python blank_contaminants_report.py      # Step 8 - re-run after blank correction
```

## Output files

```
output/
|-- peak_matrix_raw.csv                  features x samples (after feature detection)
|-- peak_matrix_blank_corrected.csv      after blank correction (full feature set)
|-- peak_matrix_processed.csv            normalized + log2 + scaled (samples x features)
|                                        used by: HCA, volcano plot
|-- peak_matrix_processed_pca.csv        same as above with prevalence + EXCLUSION_LIST applied
|                                        used by: PCA
|
|-- Provenance / audit trail
|-- feature_metadata.csv                 mean RT, m/z, compound_name (from most abundant
|                                        signal), cluster spread (rt_min/max/std,
|                                        mz_min/max/std), n_samples_detected, n_contributing_peaks
|-- feature_name_map.csv                 feature_id -> compound_name lookup table
|                                        (used by PCA / HCA / volcano when FEATURE_LABEL='name')
|-- feature_peak_log.csv                 one row per raw peak: feature_id, sample (= source file),
|                                        ref_mz, rt_raw (pre-alignment), rt_aligned, rt_shift,
|                                        area, selected (True = used in matrix; False = replaced
|                                        by higher-area duplicate from same sample in cluster)
|-- rt_alignment_shifts.csv              median RT shift applied per sample (0.0 when disabled)
|-- blank_features.csv                   max_blank_area, blank_rt, blank_mz per feature
|                                        (blank_rt/mz = RT/m/z of the best RT-matched blank peak)
|-- sample_groups.csv                    group assignments for plot colouring
|-- features_removed_blank.csv           blank filter removals: mean_rt, mean_mz,
|                                        mean_sample_area, max_blank_area, fold_change
|-- features_rescued_mz_gate.csv         features that had an RT-matched blank peak
|                                        but were retained because |Δm/z| > BLANK_MZ_TOLERANCE
|                                        (only when BLANK_USE_MZ=True and any are rescued)
|-- features_removed_exclusion.csv       exclusion list removals: mean_rt, matched_exclusion_rt,
|                                        rt_deviation  (only when EXCLUSION_LIST is non-empty)
|-- features_removed_prevalence_hca.csv  prevalence filter removals for HCA matrix
|                                        (only when MIN_PREVALENCE_HCA > 0)
|-- features_removed_prevalence_pca.csv  prevalence filter removals for PCA matrix
|                                        (only when MIN_PREVALENCE_PCA > 0)
|-- features_removed_zero_variance.csv   features dropped by scaling (zero variance after log2)
|                                        (only when any are removed)
|
|-- compound_classes.csv                 per-feature PubChem classification:
|                                        feature_id, compound_name, classification_status,
|                                        pubchem_cid, iupac_name, molecular_formula,
|                                        kingdom, superclass, class, subclass, direct_parent,
|                                        npclassifier_pathway/superclass/class
|                                        (only written when RUN_COMPOUND_CLASSIFICATION = True)
|-- feature_metadata_enriched.csv        feature_metadata.csv with class columns joined;
|                                        use for class-coloured or class-labelled plots
|-- pubchem_cache.json                   local API response cache; delete to force re-fetch
|
|-- prevalence_summary.csv               per-prevalence-level feature counts
|                                        (n_samples_detected, prevalence, n_features)
|
|-- pca_scores.csv                       sample scores for all N_COMPONENTS
|-- pca_loadings.csv                     feature loadings for all N_COMPONENTS
|-- pca_variance.csv                     explained variance per component
|-- hca_sample_order.csv                 sample dendrogram leaf order
|-- hca_feature_order.csv               feature dendrogram leaf order
|-- volcano_<A>_vs_<B>.csv              per-feature log2FC, p-value, adj. p-value
|-- top_features_analysis.csv            top N features by PCA loading magnitude,
|                                        with compound name, RT, area per sample,
|                                        and PC1/PC2(/PC3) loadings
|-- blank_contaminants_report.csv        blank-removed features with compound name,
|                                        RT, m/z, max blank area, and sample list
|-- plots/
    |-- class_pie_<column>_<group>.png   compound class pie chart per column x group
    |                                    (only when RUN_CLASS_PLOTS = True)
    |-- prevalence_histogram.png         feature prevalence distribution; one bar per k/N
    |                                    detection level; reference lines for PCA/HCA thresholds
    |-- pca_scores.png                   scores scatter plot with group ellipses
    |-- pca_loadings.png                 loadings scatter plot (top features labelled)
    |-- pca_loadings_bar.png             grouped loading bar chart (top N features)
    |-- pca_scores_3d.html               interactive 3D scores (N_COMPONENTS = 3 only)
    |-- pca_loadings_3d.html             interactive 3D loadings (N_COMPONENTS = 3 only)
    |-- hca_heatmap.png                  bidirectional clustered heatmap
    |-- volcano_<A>_vs_<B>.png           volcano plot per pairwise comparison
```

## Key parameters reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `RT_MARGIN` | `0.05` | Minutes; peaks within this window -> same feature |
| `USE_MZ` | `False` | Also cluster by m/z (for dense spectra) |
| `ALIGN_RT` | `True` | Apply median RT-shift correction before feature detection |
| `MIN_PREVALENCE_PCA` | `0.5` | Min fraction of samples a feature must be detected in (area > 0) for PCA; 0.0 = disabled |
| `MIN_PREVALENCE_HCA` | `0.0` | Same filter for HCA; disabled by default |
| `MIN_PREVALENCE_VOLCANO` | `0.0` | Same filter for volcano; keep at 0.0 — group-specific features (present in one group, absent in other) are the most relevant findings |
| `FOLD_CHANGE_THRESHOLD` | `3.0` | Minimum sample/blank ratio to retain a feature |
| `BLANK_USE_MZ` | `False` | Also require m/z proximity for a blank peak to count as a match; prevents a different compound eluting at the same RT in the blank from incorrectly removing a sample feature |
| `BLANK_MZ_TOLERANCE` | `0.005` | Da — maximum \|feature_mz − blank_mz\| accepted when `BLANK_USE_MZ = True` |
| `EXCLUSION_LIST` | `[]` | RT values (min) of known compounds to exclude from PCA only |
| `EXCLUSION_RT_MARGIN` | `0.05` | +- window for exclusion list matching |
| `NORMALIZATION` | `"sum"` | Per-sample signal correction method |
| `LOG_BASE` | `2` | Log transformation base (2 = log2, conventional in metabolomics) |
| `SCALING` | `"pareto"` | Per-feature scaling; Pareto recommended for untargeted GC-MS |
| `N_COMPONENTS` | `2` | PCA components; set to `3` for interactive 3D plots |
| `PCA_ELLIPSE` | `True` | Draw 95% confidence ellipses per group |
| `PCA_TOP_LOADINGS` | `10` | Features labelled in the loadings scatter plot |
| `PCA_BAR_TOP` | `10` | Features shown in the loading bar chart |
| `HCA_LINKAGE` | `"ward"` | Linkage method for HCA (ward requires euclidean metric) |
| `HCA_METRIC` | `"euclidean"` | Distance metric for HCA |
| `HCA_CMAP` | `"vlag"` | Colormap for heatmap (diverging, suited to scaled data) |
| `HCA_MAX_FEATURE_LABELS` | `50` | Show feature axis labels when n_features <= this |
| `HCA_CLASS_ANNOTATION_COLUMNS` | `[]` | Columns from enriched metadata to show as colored annotation strips alongside the feature dendrogram; e.g. `["superclass", "npclassifier_pathway"]` |
| `COMPOUND_NAME_COL` | `"Name"` | Column in TraceFinder CSVs holding compound names; `""` to disable |
| `FEATURE_LABEL` | `"id"` | `"id"` = use feature_id in plots; `"name"` = use compound name (falls back to feature_id) |
| `VOLCANO_COMPARISONS` | `"all"` | Group pairs to compare; "all" runs all pairwise |
| `VOLCANO_FC_THRESHOLD` | `1.0` | log2 fold-change cutoff (1.0 = 2-fold) |
| `VOLCANO_P_THRESHOLD` | `0.05` | BH-adjusted p-value significance threshold |
| `VOLCANO_TOP_LABELS` | `10` | Most significant features labelled in the volcano plot |
| `RUN_COMPOUND_CLASSIFICATION` | `True` | Run the PubChem classification step in pipeline.py; set to `False` to skip entirely |
| `RUN_CLASS_PLOTS` | `True` | Run the compound class pie chart step in pipeline.py; set to `False` to skip |
| `CLASS_PIE_COLUMNS` | `["superclass", "npclassifier_pathway"]` | Columns from enriched metadata to visualise as pie charts |
| `CLASS_PIE_GROUPS` | `"separate"` | `"separate"` = one pie per group; `"combined"` = one pie total; or a list of group names |
| `CLASS_PIE_MIN_FRACTION` | `0.02` | Slices < this fraction merged into "Other" |
| `CLASS_PIE_DETECTED_ONLY` | `True` | Count only features detected in the group; `False` = count all features |
| `CLASS_HIGHLIGHT` | `[]` | List of `{"column", "value", "color"}` dicts; matching features highlighted in PCA loadings and volcano |
| `CLASS_LABEL_COLUMN` | `""` | Column whose value is appended as `[class]` to feature labels in PCA loadings, bar chart, and volcano |
| `PUBCHEM_CACHE_ONLY` | `False` | `True` = build output from local cache only, no network requests; compounds not yet cached are skipped. `False` = fetch missing entries from PubChem as normal |
| `PUBCHEM_USER_AGENT` | `"TF_NT_pipeline/1.0 (...)"` | User-Agent header sent with all PubChem requests; replace `YOUR_EMAIL_HERE` |
| `PUBCHEM_RATE_LIMIT_DELAY` | `0.35` | Seconds between PubChem requests (~2.9/sec, under the 5/sec limit) |
| `CLASSYFIRE_RATE_LIMIT_DELAY` | `1.0` | Seconds between ClassyFire requests (no stated limit; 1.0 s is conservative) |
| `NPCLASSIFIER_RATE_LIMIT_DELAY` | `1.0` | Seconds between NPClassifier requests (no stated limit; 1.0 s is conservative) |
| `PUBCHEM_CACHE_FILE` | `"output/pubchem_cache.json"` | Local JSON cache; delete to force a full re-fetch |
