# TF_NT - Nontargeted GC-MS Metabolomics Pipeline for Tracefinder .csv exports

A modular, configurable data processing pipeline for nontargeted GC-MS metabolomics data exported from **Thermo TraceFinder** (Orbitrap MS).

## What it does

| Step | Script | Input -> Output |
|------|--------|----------------|
| 1 | `data_import.py` | TraceFinder CSVs -> `peak_matrix_raw.csv` |
| 2 | `blank_correction.py` | raw matrix -> `peak_matrix_blank_corrected.csv` |
| 3 | `normalization.py` | blank-corrected -> `peak_matrix_processed.csv` |
| 4 | `pca.py` | processed matrix -> scores/loadings plots & CSVs |

All outputs are written to the `output/` folder.

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
RT_MARGIN            = 0.05      # RT clustering window (minutes)
FOLD_CHANGE_THRESHOLD = 3.0      # blank correction stringency
NORMALIZATION        = "sum"     # "sum" | "median" | "none"
SCALING              = "pareto"  # "pareto" | "auto" | "none"
N_COMPONENTS         = 2         # PCA components (use 3 for interactive 3D plots)
```

### 3. Run the pipeline

```bash
python pipeline.py
```

Or run individual steps independently (useful when re-running after parameter changes):

```bash
python data_import.py        # Step 1 - re-run if data or RT_MARGIN changes
python blank_correction.py   # Step 2 - re-run if FOLD_CHANGE_THRESHOLD or EXCLUSION_LIST changes
python normalization.py      # Step 3 - re-run if NORMALIZATION or SCALING changes
python pca.py                # Step 4 - re-run if N_COMPONENTS or plot settings change
```

## Output files

```
output/
|---- peak_matrix_raw.csv               # features x samples (after feature detection)
|---- peak_matrix_blank_corrected.csv   # after blank correction & exclusion list
|---- peak_matrix_processed.csv         # normalized, log-transformed, scaled (samples x features)
|---- feature_metadata.csv              # mean RT and m/z per feature
|---- blank_features.csv                # max blank area per feature
|---- sample_groups.csv                 # group assignments for PCA colouring
|---- pca_scores.csv / pca_loadings.csv / pca_variance.csv
|---- features_removed_blank.csv        # audit log - blank filter removals
|---- features_removed_exclusion.csv    # audit log - exclusion list removals
|---- plots/
    |---- pca_scores.png
    |---- pca_loadings.png
    |---- pca_scores_3d.html            # interactive - only when N_COMPONENTS = 3
    |---- pca_loadings_3d.html          # interactive - only when N_COMPONENTS = 3
```

## Key parameters reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `RT_MARGIN` | `0.05` | Minutes; peaks within this window -> same feature |
| `USE_MZ` | `False` | Also cluster by m/z (for dense spectra) |
| `ALIGN_RT` | `True` | Apply median RT-shift correction before detection |
| `FOLD_CHANGE_THRESHOLD` | `3.0` | Minimum sample/blank ratio to retain a feature |
| `EXCLUSION_LIST` | `[]` | RT values (min) of known interferences to always remove |
| `EXCLUSION_RT_MARGIN` | `0.05` | +- window for exclusion list matching |
| `NORMALIZATION` | `"sum"` | Per-sample signal correction |
| `LOG_BASE` | `2` | Log transformation base |
| `SCALING` | `"pareto"` | Per-feature scaling before PCA |
| `N_COMPONENTS` | `2` | PCA components; set to `3` for interactive 3D plots |
