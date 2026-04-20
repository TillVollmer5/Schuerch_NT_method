# TF_NT - Nontargeted GC-MS Metabolomics Pipeline for TraceFinder CSV exports

A modular, configurable data processing pipeline for nontargeted GC-MS metabolomics
data exported from **Thermo TraceFinder** (Orbitrap MS).

## What it does

| Step | Script | Input → Output |
|------|--------|----------------|
| 1 | `data_import.py` | TraceFinder CSVs → `peak_matrix_raw.csv` |
| 2 | `blank_correction.py` | raw matrix → `peak_matrix_blank_corrected.csv` |
| 2b | `prevalence_histogram.py` | blank-corrected → `prevalence_histogram.png` + `prevalence_summary.csv` |
| 2c | `compound_classification.py` | `feature_metadata.csv` → `compound_classes.csv` + `feature_metadata_enriched.csv` *(optional; `RUN_COMPOUND_CLASSIFICATION`)* |
| 2d | `compound_class_plots.py` | `feature_metadata_enriched.csv` → `class_pie_*.png` *(optional; `RUN_CLASS_PLOTS`)* |
| 3 | `normalization.py` | blank-corrected → three analysis-specific matrices (see below) |
| 4 | `pca.py` | `peak_matrix_processed_pca.csv` → scores/loadings plots & CSVs |
| 5 | `hca.py` | `peak_matrix_processed_hca.csv` → clustered heatmap + dendrogram order CSVs |
| 5b | `hca_dendrogram.py` | `peak_matrix_processed_hca.csv` → `hca_dendrogram.html` *(optional; `RUN_HCA_DENDROGRAM`)* |
| 6 | `volcano.py` | `peak_matrix_processed_volcano.csv` → volcano plots & results tables |
| 7 | `top_features_analysis.py` | `pca_loadings.csv` + raw CSVs → `top_features_analysis.csv` |
| 7b | `targeted_boxplots.py` | `peak_matrix_blank_corrected.csv` → `targeted_boxplots.png` *(optional; `RUN_TARGETED_BOXPLOTS`)* |
| 8 | `blank_contaminants_report.py` | provenance files → `blank_contaminants_report.csv` |
| 9 | `classification.py` | `peak_matrix_processed_hca.csv` + metadata + references → `classification.csv` |

Each pipeline run creates a **timestamped subfolder** inside `output/`:

```
output/
└── 2026-04-14_15-30-00/      ← one folder per run
    ├── pipeline.log           ← full terminal output of that run
    ├── plots/
    ├── peak_matrix_raw.csv
    └── ...
```

The run folder name is `YYYY-MM-DD_HH-MM-SS` (local time at run start).
All output files and plots for that run are written there.
`output/pubchem_cache.json` is the only file that stays in the base `output/`
folder — it is a shared API cache reused across all runs.

### Step 3 output matrices

`normalization.py` now produces **three separate matrices**, each pre-processed
with settings appropriate for its downstream analysis:

| File | Used by | Settings |
|------|---------|---------|
| `peak_matrix_processed_pca.csv` | PCA | exclusion-list filtered; `NORMALIZATION_PCA`, `LOG_BASE_PCA`, `SCALING_PCA` |
| `peak_matrix_processed_hca.csv` | HCA + dendrogram + classification | full feature set; `NORMALIZATION_HCA`, `LOG_BASE_HCA`, `SCALING_HCA` |
| `peak_matrix_processed_volcano.csv` | Volcano plot | full feature set; `NORMALIZATION_VOLCANO`, `LOG_BASE_VOLCANO`; scaling always `"none"` |

Per-analysis overrides default to `None` (inherit global `NORMALIZATION` / `LOG_BASE` / `SCALING`).

### Note on the exclusion list and targeted list

- **`EXCLUSION_LIST`** — known compounds excluded from the *PCA* matrix only. HCA and volcano always retain them so they remain statistically testable.
- **`TARGETED_LIST`** — compounds for individual boxplot panels. Independent from `EXCLUSION_LIST`; can include any compound by RT+m/z or explicit `feature_id`.

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
EXCLUSION_LIST        = []       # compounds to withhold from PCA only
TARGETED_LIST         = []       # compounds for targeted boxplots (falls back to EXCLUSION_LIST)

# Global normalization defaults (each analysis can override via NORMALIZATION_<X> etc.)
NORMALIZATION         = "pqn"    # "pqn" | "sum" | "median" | "none"
LOG_BASE              = 2        # 2 | 10 | math.e | "sqrt" | "cbrt" | "none"
SCALING               = "pareto" # "pareto" | "auto" | "vast" | "range" | "level" | "none"

# Per-analysis overrides (None = inherit global)
NORMALIZATION_VOLCANO = None     # inherits NORMALIZATION
SCALING_VOLCANO       = "none"   # MUST stay "none" — scaling distorts fold-change
NORMALIZATION_BOXPLOT = "none"   # raw areas recommended for targeted inspection
LOG_BASE_BOXPLOT      = "none"
SCALING_BOXPLOT       = "none"

# Statistical tests
STAT_TEST_VOLCANO     = "mannwhitney"  # "mannwhitney" | "ttest" | "ttest_equal" | "kruskal"
STAT_TEST_BOXPLOT     = "mannwhitney"  # above + "wilcoxon" | "spearman" | "pearson" | "anova"

N_COMPONENTS          = 2        # PCA components (use 3 for interactive 3D plots)
VOLCANO_COMPARISONS   = "all"    # "all" for all pairwise, or [(groupA, groupB), ...]
```

### 3. Run the pipeline

```bash
python pipeline.py
```

Or run individual steps independently (useful when re-running after parameter changes):

```bash
python data_import.py              # Step 1  - re-run if data or RT_MARGIN changes
python blank_correction.py         # Step 2  - re-run if FOLD_CHANGE_THRESHOLD changes
python prevalence_histogram.py     # Step 2b - review detection rates after blank correction
python compound_classification.py  # Step 2c - refresh PubChem annotations (uses cache)
python compound_class_plots.py     # Step 2d - re-run after classification or CLASS_PIE_* changes
python normalization.py            # Step 3  - re-run if any NORMALIZATION_*/LOG_BASE_*/SCALING_* changes
python pca.py                      # Step 4  - re-run if PCA settings change
python hca.py                      # Step 5  - re-run if HCA_LINKAGE or HCA_METRIC changes
python hca_dendrogram.py           # Step 5b - interactive HTML dendrogram (optional)
python volcano.py                  # Step 6  - re-run if VOLCANO_* thresholds or STAT_TEST_VOLCANO changes
python top_features_analysis.py    # Step 7  - re-run after PCA; optional --n flag
python targeted_boxplots.py        # Step 7b - re-run if TARGETED_LIST or STAT_TEST_BOXPLOT changes
python blank_contaminants_report.py  # Step 8 - re-run after blank correction
python classification.py             # Step 9 - re-run after normalization or reference files change
```

## Output files

```
output/
|-- pubchem_cache.json               shared API cache; persists across all runs
|
└── 2026-04-14_15-30-00/             one folder per pipeline run (YYYY-MM-DD_HH-MM-SS)
    |-- pipeline.log                 full terminal output captured during this run
    |-- peak_matrix_raw.csv                  features x samples (after feature detection)
|-- peak_matrix_blank_corrected.csv      after blank correction (full feature set)
|-- peak_matrix_processed_pca.csv        normalized+log+scaled (PCA settings); prevalence
|                                        filter + EXCLUSION_LIST applied; used by: PCA
|-- peak_matrix_processed_hca.csv        normalized+log+scaled (HCA settings); full feature
|                                        set; used by: HCA, dendrogram, classification
|-- peak_matrix_processed_volcano.csv    normalized+log, no scaling; full feature set;
|                                        used by: volcano plot
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
|-- features_removed_prevalence_pca.csv  prevalence filter removals for PCA matrix
|-- features_removed_prevalence_volcano.csv  prevalence filter removals for volcano matrix
|                                        (only written when the respective MIN_PREVALENCE_* > 0)
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
|
|-- prevalence_summary.csv               per-prevalence-level feature counts
|                                        (n_samples_detected, prevalence, n_features)
|
|-- pca_scores.csv                       sample scores for all N_COMPONENTS
|-- pca_loadings.csv                     feature loadings for all N_COMPONENTS
|-- pca_variance.csv                     explained variance per component
|-- hca_sample_order.csv                 sample dendrogram leaf order
|-- hca_feature_order.csv               feature dendrogram leaf order
|-- volcano_<A>_vs_<B>.csv              per-feature log2FC, p-value, adj. p-value, direction
|-- top_features_analysis.csv            top N features by group-separation score
|                                        (between-group scatter in PC score space),
|                                        with compound name, RT, area per sample,
|                                        and PC1/PC2(/PC3) loadings
|-- targeted_boxplot_data.csv            extracted values for TARGETED_LIST compounds
|                                        after boxplot normalization; one row per compound
|-- blank_contaminants_report.csv        blank-removed features with compound name,
|                                        RT, m/z, max blank area, and sample list
|-- classification.csv                   one row per HCA feature; class (1-4),
|                                        blank-corrected area per sample, and note.
|                                        Class 1 = matched to DATA/references/references.csv
|                                        Class 2 = matched to DATA/references/farn1-46*.csv
|                                           (or SI≥800, HRF≥90, |ΔRI|≤50)
|                                        Class 3 = SI≥700, HRF≥80, |ΔRI|≤100
|                                        Class 4 = all remaining features
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
    |-- class_highlight_legend.png       standalone color key for CLASS_HIGHLIGHT rules
    |                                    (only written when CLASS_HIGHLIGHT is non-empty)
    |-- hca_heatmap.png                  bidirectional clustered heatmap
    |-- hca_class_legend.png             standalone color key for HCA annotation strips
    |                                    and sample groups
    |-- volcano_<A>_vs_<B>.png           volcano plot per pairwise comparison
    |-- targeted_boxplots.png            one panel per TARGETED_LIST compound
```

## Key parameters reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `RT_MARGIN` | `0.05` | Minutes; peaks within this window -> same feature |
| `USE_MZ` | `True` | Also cluster by m/z (for dense spectra) |
| `ALIGN_RT` | `True` | Apply median RT-shift correction before feature detection |
| `MIN_PREVALENCE_PCA` | `0.35` | Min fraction of samples a feature must be detected in (area > 0) for PCA; 0.0 = disabled |
| `MIN_PREVALENCE_HCA` | `0.45` | Same filter for HCA; set > 0 to drop sparse features from the heatmap |
| `MIN_PREVALENCE_VOLCANO` | `0.0` | Same filter for volcano; keep at 0.0 — group-specific features (present in one group, absent in other) are the most relevant findings |
| `FOLD_CHANGE_THRESHOLD` | `3.0` | Minimum sample/blank ratio to retain a feature |
| `BLANK_USE_MZ` | `True` | Also require m/z proximity for a blank peak to count as a match; prevents a different compound eluting at the same RT in the blank from incorrectly removing a sample feature |
| `BLANK_MZ_TOLERANCE` | `0.0005` | Da — maximum \|feature_mz − blank_mz\| accepted when `BLANK_USE_MZ = True` |
| `EXCLUSION_LIST` | `[]` | Compounds to withhold from PCA; matched by RT ± `EXCLUSION_RT_MARGIN` |
| `TARGETED_LIST` | `[]` | Compounds for targeted boxplots; falls back to `EXCLUSION_LIST` when empty |
| `EXCLUSION_RT_MARGIN` | `0.05` | ± RT window (min) for exclusion/targeted matching |
| `NORMALIZATION` | `"pqn"` | Global per-sample normalization: `"pqn"` \| `"sum"` \| `"median"` \| `"none"` |
| `LOG_BASE` | `2` | Global log transform: `2` \| `10` \| `math.e` \| `"sqrt"` \| `"cbrt"` \| `"none"` |
| `SCALING` | `"pareto"` | Global per-feature scaling: `"pareto"` \| `"auto"` \| `"vast"` \| `"range"` \| `"level"` \| `"none"` |
| `NORMALIZATION_PCA` / `_HCA` / `_VOLCANO` / `_BOXPLOT` | `None` | Per-analysis normalization override; `None` = inherit global |
| `LOG_BASE_PCA` / `_HCA` / `_VOLCANO` / `_BOXPLOT` | `None` | Per-analysis log override; `None` = inherit global |
| `SCALING_PCA` / `_HCA` / `_VOLCANO` / `_BOXPLOT` | `None` | Per-analysis scaling override; `None` = inherit global |
| `STAT_TEST_VOLCANO` | `"mannwhitney"` | Per-feature p-value test in volcano: `"mannwhitney"` \| `"ttest"` \| `"ttest_equal"` \| `"kruskal"` |
| `STAT_TEST_BOXPLOT` | `"mannwhitney"` | Group comparison test in boxplots: above + `"wilcoxon"` \| `"spearman"` \| `"pearson"` \| `"anova"` |
| `TARGETED_BOXPLOT_ROWS` | `3` | Number of rows in the boxplot grid |
| `RUN_TARGETED_BOXPLOTS` | `True` | Run targeted boxplot step in pipeline |
| `N_COMPONENTS` | `2` | PCA components; set to `3` for interactive 3D plots |
| `PCA_ELLIPSE` | `True` | Draw 95% confidence ellipses per group |
| `PCA_TOP_LOADINGS` | `10` | Features labelled in the loadings scatter plot |
| `PCA_BAR_TOP` | `10` | Features shown in the loading bar chart |
| `HCA_LINKAGE` | `"ward"` | Linkage method for HCA (ward requires euclidean metric) |
| `HCA_METRIC` | `"euclidean"` | Distance metric for HCA |
| `HCA_CMAP` | `"vlag"` | Colormap for heatmap (diverging, suited to scaled data) |
| `HCA_MAX_FEATURE_LABELS` | `50` | Show feature axis labels when n_features <= this |
| `HCA_CLASS_ANNOTATION_COLUMNS` | `[]` | Columns from enriched metadata to show as colored annotation strips alongside the feature dendrogram; e.g. `["superclass", "npclassifier_pathway"]` |
| `RUN_HCA_DENDROGRAM` | `True` | Generate the interactive HTML dendrogram (Step 5b); set to `False` to skip |
| `COMPOUND_NAME_COL` | `"Name"` | Column in TraceFinder CSVs holding compound names; `""` to disable |
| `FEATURE_LABEL` | `"id"` | `"id"` = use feature_id in plots; `"name"` = use compound name (falls back to feature_id) |
| `VOLCANO_COMPARISONS` | `"all"` | Group pairs to compare; `"all"` runs all pairwise |
| `VOLCANO_FC_THRESHOLD` | `1.0` | log2 fold-change cutoff (1.0 = 2-fold) |
| `VOLCANO_P_THRESHOLD` | `0.05` | BH-adjusted p-value significance threshold |
| `VOLCANO_TOP_LABELS` | `30` | Most significant features labelled in the volcano plot |
| `RUN_COMPOUND_CLASSIFICATION` | `True` | Run the PubChem classification step in pipeline.py; set to `False` to skip entirely |
| `RUN_CLASS_PLOTS` | `True` | Run the compound class pie chart step in pipeline.py; set to `False` to skip |
| `CLASS_PIE_COLUMNS` | `["superclass", "npclassifier_pathway"]` | Columns from enriched metadata to visualise as pie charts |
| `CLASS_PIE_GROUPS` | `"separate"` | `"separate"` = one pie per group; `"combined"` = one pie total; or a list of group names |
| `CLASS_PIE_MIN_FRACTION` | `0.02` | Slices < this fraction merged into "Other" |
| `CLASS_PIE_DETECTED_ONLY` | `True` | Count only features detected in the group; `False` = count all features |
| `CLASS_HIGHLIGHT` | `[]` | Ordered list of `{"column", "value", "color"}` dicts. Rules are applied in order; **later rules overwrite earlier ones**, so put broad classes first and specific subclasses last. Matching features appear as colored dots in PCA loadings/bar chart, colored ring outlines on significant volcano dots, and colored tick labels in the HCA bar chart. The same hex colors are reused in the HCA annotation strip legend so all plots are consistent. A standalone `class_highlight_legend.png` is written automatically. |
| `CLASS_LABEL_COLUMN` | `""` | Column whose value is appended as `[class]` to feature labels in PCA loadings, bar chart, and volcano |
| `PUBCHEM_CACHE_ONLY` | `False` | `True` = build output from local cache only, no network requests; compounds not yet cached are skipped. `False` = fetch missing entries from PubChem as normal |
| `PUBCHEM_USER_AGENT` | `"TF_NT_pipeline/1.0 (...)"` | User-Agent header sent with all PubChem requests; replace `YOUR_EMAIL_HERE` |
| `PUBCHEM_RATE_LIMIT_DELAY` | `1` | Seconds between PubChem requests (conservative; increase on frequent 503 responses) |
| `CLASSYFIRE_RATE_LIMIT_DELAY` | `3` | Seconds between ClassyFire requests (no stated limit; 3 s avoids HTTP 429) |
| `NPCLASSIFIER_RATE_LIMIT_DELAY` | `3` | Seconds between NPClassifier requests (no stated limit; 3 s is conservative) |
| `PUBCHEM_CACHE_FILE` | `"output/pubchem_cache.json"` | Local JSON cache shared across all runs; stays in base `output/`, not the timestamped run folder; delete to force a full re-fetch |
