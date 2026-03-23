# TF_NT - Nontargeted GC-MS Metabolomics Pipeline for TraceFinder CSV exports

A modular, configurable data processing pipeline for nontargeted GC-MS metabolomics
data exported from **Thermo TraceFinder** (Orbitrap MS).

## What it does

| Step | Script | Input -> Output |
|------|--------|----------------|
| 1 | `data_import.py` | TraceFinder CSVs -> `peak_matrix_raw.csv` |
| 2 | `blank_correction.py` | raw matrix -> `peak_matrix_blank_corrected.csv` |
| 3 | `normalization.py` | blank-corrected -> `peak_matrix_processed.csv` + `peak_matrix_processed_pca.csv` |
| 4 | `pca.py` | `peak_matrix_processed_pca.csv` -> scores/loadings plots & CSVs |
| 5 | `hca.py` | `peak_matrix_processed.csv` -> clustered heatmap + dendrogram order CSVs |
| 6 | `volcano.py` | `peak_matrix_blank_corrected.csv` -> volcano plots & results tables |

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
python data_import.py        # Step 1 - re-run if data or RT_MARGIN changes
python blank_correction.py   # Step 2 - re-run if FOLD_CHANGE_THRESHOLD changes
python normalization.py      # Step 3 - re-run if NORMALIZATION, SCALING, or EXCLUSION_LIST changes
python pca.py                # Step 4 - re-run if PCA settings change
python hca.py                # Step 5 - re-run if HCA_LINKAGE or HCA_METRIC changes
python volcano.py            # Step 6 - re-run if VOLCANO_* thresholds change
```

## Output files

```
output/
|-- peak_matrix_raw.csv                  features x samples (after feature detection)
|-- peak_matrix_blank_corrected.csv      after blank correction (full feature set)
|-- peak_matrix_processed.csv            normalized + log2 + scaled (samples x features)
|                                        used by: HCA, volcano plot
|-- peak_matrix_processed_pca.csv        same as above with EXCLUSION_LIST applied
|                                        used by: PCA
|-- feature_metadata.csv                 mean RT and m/z per feature
|-- blank_features.csv                   max blank area per feature
|-- sample_groups.csv                    group assignments for plot colouring
|-- features_removed_blank.csv           audit log - blank filter removals
|-- features_removed_exclusion.csv       audit log - exclusion list removals (if any)
|-- pca_scores.csv                       sample scores for all N_COMPONENTS
|-- pca_loadings.csv                     feature loadings for all N_COMPONENTS
|-- pca_variance.csv                     explained variance per component
|-- hca_sample_order.csv                 sample dendrogram leaf order
|-- hca_feature_order.csv               feature dendrogram leaf order
|-- volcano_<A>_vs_<B>.csv              per-feature log2FC, p-value, adj. p-value
|-- plots/
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
| `FOLD_CHANGE_THRESHOLD` | `3.0` | Minimum sample/blank ratio to retain a feature |
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
| `VOLCANO_COMPARISONS` | `"all"` | Group pairs to compare; "all" runs all pairwise |
| `VOLCANO_FC_THRESHOLD` | `1.0` | log2 fold-change cutoff (1.0 = 2-fold) |
| `VOLCANO_P_THRESHOLD` | `0.05` | BH-adjusted p-value significance threshold |
| `VOLCANO_TOP_LABELS` | `10` | Most significant features labelled in the volcano plot |
